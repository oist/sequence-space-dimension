#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Libraries
import os
import os.path
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
import pandas as pd

## List of genomes
eggNOG_df = pd.read_csv('./1236_members.tsv', sep='\t')
eggNOG_df.columns =['Tax_id', 'Ortho_id', 'Number_of_genomes', 'Number_of_genes', 'List_of_genes', 'List_of_tax']
eggNOG_df

ORG=[]
GN=set()
for index, row in eggNOG_df.iterrows():
    for q in row['List_of_genes'].split(','):
        GN.append(q)
    for g in row['List_of_tax'].split(','):
        ORG.add(int(g))
        
print(ORG[0])
refseq_data.loc[refseq_data['taxid'] == ORG[0]]
genome_data = pd.read_csv('tax2fname.txt', sep='\t')
genome_data.columns =['Tax_id', 'Assembly_id', 'Organism']
genome_data = genome_data.loc[genome_data['Tax_id'].isin(ORG)]
assemblies = set(genome_data['Assembly_id'])

names = []
for assembly in assemblies:
    names.append(assembly.split('_')[0]+'_'+assembly.split('_')[1])

genome_links = open('Genome_links_filtered.txt','w')

refseq_data = pd.read_csv('assembly_summary_refseq.txt', sep='\t')

for index, row in refseq_data.iterrows():
    if row['assembly_accession'] in names:
        if row['taxid'] in ORG:
            link = row['ftp_path']+'/'+row['ftp_path'].split('/')[9]+'_genomic.gbff.gz'
            genome_links.write(link+'\n')

genome_links.close()
## Data from eggNOG
RefSeq_df = pd.read_csv('../All/assembly_summary_refseq.txt', sep='\t') #genes info

## Sequences extraction
#downloaded genome parsing with genes


eggNOG_df = pd.read_csv('./1236_members.tsv', sep='\t')
eggNOG_df.columns =['Tax_id', 'Ortho_id', 'Number_of_genomes', 'Number_of_genes', 'List_of_genes', 'List_of_tax']
eggNOG_df

ORG=[]
GN=[]
for index, row in eggNOG_df.iterrows():
    for q in row['List_of_genes'].split(','):
        GN.append(q)
    for g in row['List_of_tax'].split(','):
        ORG.append(g)
        
tax = []
loci = []
for gene in GN:
    tax.append(gene.split('.')[0])
    loci.append(gene.split('.')[1])
    
print(len(loci))

GenesTable = open('./GenesTable_locus.tsv','w') #genes info
GenesTable.write('Assembly_id'+'\t'+'Gene_tag'+'\t'+'Protein_name'+'\t'+'translation'+'\t'+'transcript'+'\n')

Statistics = open('Genome_statistics_locus.csv','w')
Statistics.write('Assembly'+'\t'+'number of eggNOGs'+'\n')

for file in os.listdir("./genomes"):
    Statistics.write(str(file))
    counter_genes = 0
    for rec in SeqIO.parse("./genomes/"+file, "genbank"):
        assembly = str(file)[:-13].split('_')[0]+'_'+str(file)[:-13].split('_')[1]
        for feature in rec.features:
            if feature.type=='CDS':
                if 'translation' in feature.qualifiers:
                    if 'old_locus_tag' not in feature.qualifiers:
                        if feature.qualifiers["locus_tag"][0] in loci:
                            counter_genes+=1
                            if feature.location.strand == 1:
                                transcript = rec.seq[feature.location.start:feature.location.end]
                            else:
                                transcript = rec.seq[feature.location.start:feature.location.end].reverse_complement()
                            GenesTable.write(assembly+'\t'+feature.qualifiers["locus_tag"][0]+'\t'+feature.qualifiers["protein_id"][0]+'\t'+str(feature.qualifiers["translation"][0])+'\t'+str(transcript)+'\n')
                    else:
                        if feature.qualifiers["old_locus_tag"][0] in loci:
                            counter_genes+=1
                            if feature.location.strand == 1:
                                transcript = rec.seq[feature.location.start:feature.location.end]
                            else:
                                transcript = rec.seq[feature.location.start:feature.location.end].reverse_complement()
                            GenesTable.write(assembly+'\t'+feature.qualifiers["old_locus_tag"][0]+'\t'+feature.qualifiers["protein_id"][0]+'\t'+str(feature.qualifiers["translation"][0])+'\t'+str(transcript)+'\n')
                        
                        
    Statistics.write('\t'+str(counter_genes)+'\n')
    print(file)
GenesTable.close()
Statistics.close()

# eggNOG clustering
genes_df = pd.read_csv('./GenesTable_locus.tsv', sep='\t')

print(genes_df.shape)
print(genes_df['Gene_tag'].unique)

gene_list = list(genes_df['Gene_tag'].unique())

aa_dict = dict(zip(genes_df['Gene_tag'],genes_df['translation']))
print(len(aa_dict))

nt_dict = dict(zip(genes_df['Gene_tag'],genes_df['transcript']))
print(len(nt_dict))

genes_df
aa_dict['GNIT_0053']
eggNOG_df = pd.read_csv('./1236_members.tsv', sep='\t')
eggNOG_df.columns =['Tax_id', 'Ortho_id', 'Number_of_genes', 'Number_of_genomes', 'List_of_genes', 'List_of_tax']
eggNOG_df        
eggNOG_df.loc[eggNOG_df['Ortho_id'] == '1RM7I']
aa_f = './eggNOGs_aa/'
nt_f = './eggNOGs_nt/'

#Statistics = open('eggNOG_statistics_locus.csv','a')
#Statistics.write('name'+'\t'+'number of eggNOGs'+'\n')

for index, row in eggNOG_df.iterrows():
    if os.path.exists(nt_f+row['Ortho_id']+'.fna') != True:
        ORG=[]
        GN=[]
        eg_counter = 0
        eggNOG_aa = open(aa_f+row['Ortho_id']+'.faa', 'w')
        eggNOG_nt = open(nt_f+row['Ortho_id']+'.fna', 'w')
        for q in row['List_of_genes'].split(','):
            GN.append(q)
        for g in row['List_of_tax'].split(','):
            ORG.append(g)
        for i in range(len(GN)):
            ob = GN[i].split('.')[1]
            if ob in gene_list:
                eggNOG_aa.write('>'+GN[i]+'\n'+str(aa_dict[ob]) +'\n')
                eggNOG_nt.write('>'+GN[i]+'\n'+str(nt_dict[ob]) +'\n')
                eg_counter += 1
        eggNOG_aa.close()
        eggNOG_nt.close()
#    Statistics.write(row['Ortho_id']+'\t'+str(eg_counter)+'\n')
#Statistics.close()
        
Statistics = open('eggNOG_stats.tsv','w')

for file in os.listdir("./eggNOGs_nt/"):
    counter_genes = 0
    Statistics.write(str(file)+'\t')
    for rec in SeqIO.parse("./eggNOGs_nt/"+file, "fasta"):
        counter_genes += 1 
    Statistics.write(str(counter_genes)+'\n')
    if counter_genes < 10:
        os.remove("./eggNOGs_aa/"+file)
        os.remove("./eggNOGs_nt/"+file)
# eggNOG alignments
# aa alignments to nt alignnments

#aa_f = './eggNOGs_aa_aligned_test/'
aa_f = './1236_raw_algs/1236/'
nt_f = './eggNOGs_nt/'
out_f = './eggNOGs_nt_aligned/'

nt_egg_dict = dict()
  
for al in os.listdir(aa_f):
    egg_name = al.split('.')[0]
    if os.path.exists(nt_f+egg_name+'.fna') == True:
        out = open(out_f+egg_name+'.fasta', 'w')

        for rec_nt in SeqIO.parse(nt_f+egg_name+'.fna', "fasta"):
            nt_egg_dict[rec_nt.id] = rec_nt.seq
    
        for al_aa in SeqIO.parse(aa_f+al, "fasta"):

            gene = al_aa.id
            if gene in nt_egg_dict.keys():
                al_nt = ''
                q = 0
                for i in range(len(al_aa.seq)):
                    if al_aa.seq[i] == '-':
                        al_nt = al_nt+'---'
                    else:
                        al_nt = al_nt+nt_egg_dict[gene][q*3:q*3+3]
                        q = q+1
                if len(al_aa.seq)*3 == len(al_nt):
                    out.write('>'+gene+'\n'+str(al_nt)+'\n')
#                else:
#                    print(al_aa.seq)
#                    print(al_nt)
#    else:
#        print(nt_f+egg_name+'.fna')
    out.close()