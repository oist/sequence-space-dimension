import argparse
import hashlib
import importlib
import io
import os
import pickle
import runpy
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "src" / "simulate_evol"
RANDOM = "simulate_random_seqs_from_matrix_dim_for_diff_nseq.py"
TREE = "simulate_tree_seqs_dim_for_diff_nseq.py"
SUBSAMPLE = "subsample_tree_seqs_dim_for_diff_nseq.py"
FLAGS = (
    (RANDOM, ("-n", "7"), "-nm", "--newmatr", True),
    (TREE, ("-n", "7"), "-nm", "--newmatr", True),
    (TREE, ("-n", "7"), "-s", "--starlike", False),
    (SUBSAMPLE, ("-f", "unused"), "-m", "--matr_to_file", False),
)


def run_script(script, args, simulation=None, sequence_space=None):
    # All CLI and file-writing code is real. Simulation dependencies are test doubles.
    # Load these once outside patch.dict, which restores the import cache afterward.
    importlib.import_module("pandas")
    importlib.import_module("tqdm")
    modules = {
        "simulate_prot_evol_lib": simulation or SimpleNamespace(),
        "sequence_space_lib": sequence_space or SimpleNamespace(),
    }
    with patch.dict(sys.modules, modules), patch.object(
        sys, "argv", [script, *args],
    ), patch.object(sys, "path", sys.path.copy()), redirect_stdout(io.StringIO()):
        runpy.run_path(str(SCRIPTS / script), run_name="__main__")


class ParsedArguments(Exception):
    pass


def parse_arguments(script, args):
    parse_args = argparse.ArgumentParser.parse_args

    def capture(parser):
        raise ParsedArguments(parse_args(parser))

    with patch.object(argparse.ArgumentParser, "parse_args", capture):
        try:
            run_script(script, args)
        except ParsedArguments as result:
            return vars(result.args[0])
    raise AssertionError("Script did not parse its arguments")


def matrix_hash(matrix):
    return hashlib.sha256(matrix.tobytes()).hexdigest()


def simulation_fixture():
    matrices_seen = []
    replacements = [np.full((2, 2), value) for value in (2, 3, 4)]
    simulation = SimpleNamespace(
        generate_matrix_AA=Mock(side_effect=replacements),
        generate_nt_seq=Mock(side_effect=lambda *args, matrix, **kw: matrices_seen.append(matrix)),
        simulate_tree_with_depth_scaled_branch_lengths=Mock(return_value=object()),
        evolve_tree=Mock(return_value=([], 0, 0, 0, None, [])),
        get_dnds_dist_dim_stats_from_seq=Mock(return_value=([0.1, 0.2], 0.15, (1, 1, 0))),
    )
    return simulation, matrices_seen, replacements


class SimulationBooleanTests(unittest.TestCase):
    def test_accepts_boolean_spellings_with_both_option_names(self):
        for script, required, short, long, _ in FLAGS:
            for wanted, spellings in (
                (True, ("True", "true", "TRUE", "T", "yes", "Y", "1")),
                (False, ("False", "false", "FALSE", "F", "no", "N", "0")),
            ):
                for option in (short, long):
                    for value in spellings:
                        with self.subTest(script=script, option=option, value=value):
                            actual = parse_arguments(script, [*required, option, value])
                            self.assertIs(actual[long[2:]], wanted)

    def test_rejects_invalid_values(self):
        for script, required, short, long, _ in FLAGS:
            for value in ("maybe", "offish", "falsee", "none", ""):
                with self.subTest(script=script, option=short, value=value):
                    messages = io.StringIO()
                    with redirect_stderr(messages), self.assertRaises(SystemExit) as error:
                        parse_arguments(script, [*required, short, value])
                    self.assertEqual(error.exception.code, 2)
                    self.assertIn(long, messages.getvalue())
                    self.assertIn("expected a boolean value", messages.getvalue())

    def test_defaults_and_numeric_arguments(self):
        for script, required, _, long, default in FLAGS:
            with self.subTest(script=script, option=long):
                counts = ["-n", "2"] if script == SUBSAMPLE else ["-r", "2"]
                parsed = parse_arguments(script, [*required, *counts, "-w", "0.25"])
                self.assertIs(parsed[long[2:]], default)
                self.assertEqual(parsed["win_size"], 0.25)
                self.assertEqual(parsed["nsamples" if script == SUBSAMPLE else "repl"], 2)
                if script != SUBSAMPLE:
                    self.assertEqual(parsed["nseq"], 7)

    def check_matrix_mode(self, regenerate):
        for script in (RANDOM, TREE):
            with self.subTest(script=script), tempfile.TemporaryDirectory() as directory:
                source = np.arange(4).reshape(2, 2)
                matrix_path = Path(directory) / "supplied.pickle"
                original = pickle.dumps([(0.5, source)])
                matrix_path.write_bytes(original)
                simulation, seen, replacements = simulation_fixture()
                run_script(script, [
                    "-n", "7", "-r", "3", "-m", str(matrix_path),
                    "-of", directory + os.sep, "-nm", str(regenerate),
                ], simulation)
                expected = replacements if regenerate else [source] * 3
                self.assertEqual([matrix_hash(m) for m in seen], [matrix_hash(m) for m in expected])
                self.assertEqual(len({id(matrix) for matrix in seen}), 3 if regenerate else 1)
                self.assertEqual(simulation.generate_matrix_AA.call_count, 3 if regenerate else 0)
                if script == TREE:
                    self.assertEqual(
                        [id(call.args[1]) for call in simulation.evolve_tree.call_args_list],
                        [id(matrix) for matrix in seen],
                    )
                self.assertEqual(matrix_path.read_bytes(), original)

    def test_false_reuses_supplied_matrix(self):
        self.check_matrix_mode(False)

    def test_true_regenerates_matrix_for_each_replicate(self):
        self.check_matrix_mode(True)

    def test_starlike_selects_tree_branch(self):
        for option, starlike, scale in (([], False, 1), (["-s", "True"], True, 0.7),
                                         (["-s", "False"], False, 1)):
            with self.subTest(option=option), tempfile.TemporaryDirectory() as directory:
                matrix_path = Path(directory) / "supplied.pickle"
                matrix_path.write_bytes(pickle.dumps([(0.5, np.eye(2))]))
                simulation, _, _ = simulation_fixture()
                run_script(TREE, [
                    "-n", "7", "-r", "1", "-m", str(matrix_path), "-of", directory + os.sep,
                    *option,
                ], simulation)
                simulation.simulate_tree_with_depth_scaled_branch_lengths.assert_called_once_with(
                    7, starlike=starlike, length_scale=scale,
                )

    def test_subsampling_matrix_output_flag(self):
        for option, writes in (([], False), (["-m", "True"], True), (["-m", "False"], False)):
            with self.subTest(option=option), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                file_list = root / "matrices.txt"
                file_list.write_text("simu_tree_f0.5_nseq10000_rep0.dist.pdistm\n")
                seqspace = SimpleNamespace(
                    submatrices=Mock(return_value=[np.eye(3)]),
                    matrix_to_array=Mock(return_value=([0.1, 0.2, 0.3], None)),
                    max_k=Mock(return_value=(1, 1, 0)),
                )
                run_script(SUBSAMPLE, [
                    "-f", str(file_list), "-s", "[3]", "-n", "1",
                    "-of", directory + os.sep, *option,
                ], sequence_space=seqspace)
                outputs = list(root.glob("*.pdistm"))
                self.assertEqual(len(outputs), int(writes))
                if writes:
                    np.testing.assert_array_equal(np.loadtxt(outputs[0], delimiter=","), np.eye(3))
                seqspace.max_k.assert_called_once()


if __name__ == "__main__":
    unittest.main()
