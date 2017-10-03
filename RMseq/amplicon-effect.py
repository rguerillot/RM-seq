#!/usr/bin/env python3

'''
    Uses python3.
    email: romain.guerillot@hotmail.fr
    Authors: Romain Guerillot, Torsten Seemann, Mark B. Schultz
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from plumbum import local
from plumbum.cmd import getorf, diffseq, grep, rm
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import os
import argparse
import tempfile
import subprocess
# import progressbar

#argument parser
parser = argparse.ArgumentParser("Translate nucleotide multi-fasta file to amino acid and compare each peptide with reference protein and output mutation annotation in a tabular format")
parser.add_argument("mfasta", type=str, help="multi-fasta file")
parser.add_argument("ref", type=str, help="reference protein (fasta)")
parser.add_argument("outdir", help="output directory")
parser.add_argument("-n", "--name", default = "sample", help="sample name to add in table (default: sample)")
parser.add_argument("-f", "--filter", default = "200", help="Minimum size of ORF in bp to annotate (default: 200)")
parser.add_argument("-w", "--wordsize", default = "5", help="The similar regions between the two sequences are found by creating a hash table of 'wordsize'd subsequences. 10 is a reasonable default. Making this value larger (20) may speed up the program slightly, but will mean that any two differences within 'wordsize' of each other will be grouped as a single region of difference. This value may be made smaller (4) to improve the resolution of nearby differences, but the program will go much slower. (default: 5)")

args = parser.parse_args()
conseq_file = os.path.abspath(args.mfasta)
reference_file = os.path.abspath(args.ref)
sample_name = args.name
outfolder = os.path.abspath(args.outdir)
ORFminsize = args.filter
wordsize = args.wordsize
ref = SeqIO.read(reference_file, "fasta")

def deleteContent(pfile):
    pfile.seek(0)
    pfile.truncate()

def get_aa_from_position(sstart, send, ref):
    if send == sstart:
        return str(ref.seq[int(sstart)-1])
    else:
        return str(ref.seq[int(sstart)-1:int(send)])
        

def reformat_diffseq_output(diffseq_output, barcode, dict_bar_nuc, dict_bar_orf):
    if 'No hits in output report file' in diffseq_output :
        clean_lines = sample_name + "\t" + "WT" + "\t" + "None" + "\t" + "None" + "\t" +  dict_bar_nuc.get(str(barcode)) + "\t" + dict_bar_orf.get(str(barcode)) +  "\n"
        raw_lines = sample_name + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + dict_bar_nuc.get(str(barcode)) + "\t" + dict_bar_orf.get(str(barcode)) +  "\n"
        return(raw_lines, clean_lines)
    out_header ="Compare and report features of two similar sequences\n" + "SeqName" + "\t" + "Start" + "\t" + "End" + "\t" + "Score" + "\t" + "Strand" + "\t" + "start" + "\t" + "end" + "\t" + "length" + "\t" + "name" + "\t" + "sequence" + "\t" + "first_feature" + "\t" + "second_feature" + "\n"
    lines = diffseq_output.replace(out_header, "") #remove header from output
    lines = lines[0:-2] # remove last empty newline
    lines = lines.split('\n')
    raw_lines = ""
    clean_lines = ""
    my_mut_nb = 0
    mutation = []
    while my_mut_nb < len(lines) :
        line1 = lines[my_mut_nb].split("\t")
        start = line1[1]
        end = line1[2]
        refaa = get_aa_from_position(start, end, ref)
        sample = line1[8].split("---")[0]
        conseqtag = line1[8].split("---")[1].split("_")[0]
        mutation.append(refaa + start + line1[9])
        my_mut_nb += 1
    raw_lines += sample +"\t" + "\t".join(line1[1::]) + "\t" +  dict_bar_nuc.get(str(barcode)) + "\t" + dict_bar_orf.get(str(barcode)) +  "\n"
    clean_lines += sample + "\t" + ", ".join(mutation) + "\t" + start + "\t" + end + "\t" +  dict_bar_nuc.get(str(barcode)) + "\t" + dict_bar_orf.get(str(barcode)) +  "\n"
    return(raw_lines, clean_lines)

def get_barcode_cluster():
    with open(outfolder + "/amplicons.cdhit.bak.clstr") as f:
        lines = f.readlines()
        bar_clus = {}
        for line in lines:
            clusternb = line.split("\t")[0]
            barcode = line.split(">")[1].split("...")[0]
            bar_clus[barcode] = clusternb
    return bar_clus

def get_barcode_fasta(my_fasta):
    fasta_sequences = SeqIO.parse(open(my_fasta),'fasta')
    bar_fasta = {}
    for fasta in fasta_sequences:
        clean_fasta_id = fasta.id.split("_")[0] # to remove _{d}
        bar_fasta[clean_fasta_id] = str(fasta.seq)
    return bar_fasta

def get_cluster_annot(barcode_cluster_dict, dict_bar_nuc, dict_bar_orf):
    fasta_sequences = SeqIO.parse(open(outfolder + "/amplicons.orf"),'fasta')
    nb_fa_records = len(list(SeqIO.parse(open(outfolder + "/amplicons.orf"),'fasta')))
    print("running: diffseq")
    print("processing " + str(nb_fa_records)+ " unique orf")
    counter = 0
    cluster_number = 0
    dict_cluster_annot = {}
    with tempfile.NamedTemporaryFile(mode = "w") as my_fasta:
        for fasta in fasta_sequences:
            fabarcode = fasta.id
            fabarcode_clean = fabarcode.split("_")[0] # to remove _{d}
            fasta.id = sample_name + "---" + fasta.id
            counter += 1
            SeqIO.write(fasta, my_fasta, "fasta")
            my_fasta.flush()
            output = subprocess.check_output(['diffseq', '-asequence'] +
                                             [str(reference_file)] +
                                             ['-bsequence'] +
                                             [str(my_fasta.name)] +
                                             ['-aoutfeat',
                                              'seq1out',
                                              '-boutfeat',
                                              'seq2out',
                                              '-wordsize',
                                              wordsize,
                                              '-rdesshow3',
                                              '-rformat',
                                              'excel',
                                              'stdout'],
                                             stderr=subprocess.STDOUT) \
                                             .decode('UTF-8')
            deleteContent(my_fasta)
            key = int(barcode_cluster_dict.get(fabarcode_clean))
            dict_cluster_annot[key] = reformat_diffseq_output(output,
                                                              fabarcode_clean,
                                                              dict_bar_nuc,
                                                              dict_bar_orf)
            cluster_number += 1
#                 bar.update(counter)
#         bar.finish()
    return dict_cluster_annot
        
def main():
    print("multi-fasta file: " + conseq_file)
    print("reference file: " + reference_file)
    print("output directory: " + outfolder)
    print("sample name: " + sample_name)
    print("keeping ORF: >= " + ORFminsize + "bp")

    # create outputdir if it doesn't exit
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    # run cd-hit-est to cluster identical amplicons
    print("running cd-hit to cluster identical amplicons")
    cdhit = subprocess.call(["cd-hit-est", "-bak", "1", "-c", "1", "-G", "0", "-aL", "1", "-i", conseq_file, "-o", outfolder + "/amplicons.cdhit"])

    # run EMBOSS getORF to find the largest ORF of each conseq
    get_ORF = getorf["-sequence", outfolder + "/amplicons.cdhit", "-outseq", outfolder + "/amplicons.orf", "-table", "11", "-minsize", ORFminsize, "-reverse", "FALSE"]
    print("running: ", get_ORF)
    get_ORF()
        
    # get dictionnary of cluster, orf and nucleotide
    barcode_cluster = get_barcode_cluster()
    barcode_orf = get_barcode_fasta(outfolder + "/amplicons.orf")
    barcode_nuc = get_barcode_fasta(conseq_file)

    # run EMBOSS diffseq on each conseq.orf entry
    try:
        rm[outfolder + "/amplicons_raw.effect"]()
        rm[outfolder + "/amplicons.effect"]()
    except:
        pass
    clust_annot = get_cluster_annot(barcode_cluster, barcode_orf, barcode_nuc)

    # write annotation files
    print("writing annotation files")
    out_raw_table = open(outfolder+"/amplicons_raw.effect", "a") #warning appending to file'
    out_raw_table.write("barcode" + "\t" + "sample" + "\t" +  "Start" + "\t" + "End" + "\t" + "Score" + "\t" + "Strand" + "\t" + "start" + "\t" + "end" + "\t" + "length" + "\t" + "name" + "\t" + "sequence" + "\t" + "first_feature" + "\t" + "second_feature" + "\t" + "orf" + "\t" + "dna" + "\n")
    out_table = open(outfolder+"/amplicons.effect", "a") #warning appending to file
    out_table.write("barcode" + "\t" + "sample" + "\t" + "aa_mutation" + "\t" + "start" + "\t" + " end" + "\t" + "orf" + "\t" + "dna" + "\n")

    print("# consensus amplicons = " + str(len(barcode_nuc)))
    print("# non-identical consensus amplicons = " + str(len(set(barcode_cluster.values()))))
    print("# non-identical consensus amplicons translated = " + str(len(barcode_orf)))
    print("# non-identical orf annotated = " + str(len(clust_annot)))


    for barcode in barcode_nuc:
        try:
            raw_annot = str(barcode) + "\t" +  clust_annot.get(int(barcode_cluster.get(str(barcode))))[0]
            annot = str(barcode) + "\t" + clust_annot.get(int(barcode_cluster.get(str(barcode))))[1]
        except TypeError: # manage .get error when there is no orf corresponding to barcode (no orf found)
            raw_annot = str(barcode) + "\t" + sample_name + "\t" + str(barcode) + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + str(barcode) + "\t" + "None" + "\t" + "None" + "\t" + "None" + "\t" + barcode_nuc.get(str(barcode)) + "\t" + "no orf found" +  "\n"
            annot = str(barcode) + "\t" + sample_name + "\t" + str(barcode) + "\t" + "NA" + "\t" + "None" + "\t" + "None" + "\t" +  barcode_nuc.get(str(barcode)) + "\t" + "no orf found" +  "\n"
        out_raw_table.write(raw_annot)
        out_table.write(annot)
    out_raw_table.close()
    out_table.close()
    os.remove("seq1out")
    os.remove("seq2out")

if __name__ == '__main__':
    main()
