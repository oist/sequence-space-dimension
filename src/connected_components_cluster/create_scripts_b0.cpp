#include <iostream>
#include <string>
#include <fstream>
using namespace std;

ifstream SIZE;
ofstream OUTPUT;
ofstream COMMANDS_DRAW;

int main(){
    SIZE.open("size_simulate_evol.txt");
    COMMANDS_DRAW.open("drawing_curve_b0_simulate_evol.txt");

    string number;
    while(SIZE >> number){
        int n;
        SIZE >> n;
        string SCRIPT, HEAD, DRAWING, CURVE;
        SCRIPT = "Gnuplot_script_b0_simulate_evol/" + number + "_script.p";
        OUTPUT.open(SCRIPT.c_str());
        HEAD = "Betti 0 for " + number + "   n = ";
        DRAWING = "Drawing_curve_b0_simulate_evol/" + number + "_drawing.png";
        CURVE = "Result_b0curve_simulate_evol/" + number + ".b0curve";
        
        OUTPUT<<"# Set the output to a png file"<<endl;
        OUTPUT<<"set terminal png"<<endl;
        OUTPUT<<"set xrange [0:1]"<<endl;
        OUTPUT<<"set yrange [0:600]"<<endl;
        OUTPUT<<"set nokey"<<endl;
        OUTPUT<<"# The file we'll write to"<<endl;
        OUTPUT<<"set output '"<< DRAWING <<"'"<<endl;
        OUTPUT<<"# The graphic title"<<endl;
        OUTPUT<<"set title '"<< HEAD << n <<"'"<<endl;
        OUTPUT<<"#plot the graphic"<<endl;
        OUTPUT<<"plot \""<<CURVE<<"\""<<' '<<"with steps lw 2";
        OUTPUT.close();
        
        COMMANDS_DRAW << "gnuplot Gnuplot_script_b0_simulate_evol/" << number << "_script.p" << endl;

    }
    SIZE.close();
    COMMANDS_DRAW.close();
    return 0;
}
