#!/home/linuxbrew/.linuxbrew/opt/python/bin/python3.7

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

#argument parser
parser = argparse.ArgumentParser("Translate nucleotide multi-fasta file to amino acid and compare each peptide with reference protein and output mutation annotation in a tabular format")
parser.add_argument("mfasta", type=str, help="Multi-fasta file")
parser.add_argument("ref", type=str, help="Nucleotide sequence of the reference gene (fasta)")
parser.add_argument("outdir", help="Output directory")
parser.add_argument("-n", "--name", default = "sample", help="Sample name to add in table (default: sample)")
parser.add_argument("-t", "--translation", default = "getorf", help="By default the pipeline will use getorf to translate amplicons to protein sequence (set the filter option). Use this option if you want to manually set the reading frame 1, 2 or 3 for translating your consensus sequence")
parser.add_argument("-f", "--filter", default = "200", help="Minimum size of ORF in bp to annotate (default: 200)")
parser.add_argument("-w", "--wordsize", default = "5", help="The similar regions between the two sequences are found by creating a hash table of 'wordsize'd subsequences. 5 is a reasonable default. Making this value larger (20) may speed up the program slightly, but will mean that any two differences within 'wordsize' of each other will be grouped as a single region of difference. This value may be made smaller (4) to improve the resolution of nearby differences, but the program will go much slower. (default: 5)")

# argument and variables
args = parser.parse_args()
conseq_file = os.path.abspath(args.mfasta)
sample_name = args.name
outfolder = os.path.abspath(args.outdir)
ORFminsize = args.filter
wordsize = args.wordsize
readingframe = args.translation


def deleteContent(pfile):
    pfile.seek(0)
    pfile.truncate()

def get_aanuc_from_position(sstart, send, ref_file):
    ref = SeqIO.read(os.path.abspath(ref_file), "fasta")
    if send == sstart:
        return str(ref.seq[int(sstart)-1])
    else:
        return str(ref.seq[int(sstart)-1:int(send)])
        
def parse_diffseq_output(diffseq_output, ref_file):
    if 'No hits in output report file' in diffseq_output:
        return(["NA", "NA", "0", "0","0"])
    out_header ="Compare and report features of two similar sequences\n" + "SeqName" + "\t" + "Start" + "\t" + "End" + "\t" + "Score" + "\t" + "Strand" + "\t" + "start" + "\t" + "end" + "\t" + "length" + "\t" + "name" + "\t" + "sequence" + "\t" + "first_feature" + "\t" + "second_feature" + "\n"
    lines = diffseq_output.replace(out_header, "") #remove header from output
    lines = lines[0:-2] # remove last empty newline
    lines = lines.split('\n')
    raw_lines = ""
    clean_lines = ""
    mut_nb = 0
    mutations = []
    starts = []
    ends = []
    while mut_nb < len(lines) :
        line1 = lines[mut_nb].split("\t")
        start = line1[1]
        starts.append(start)
        end = line1[2]
        ends.append(end)
        refaa = get_aanuc_from_position(start, end, ref_file)
        conseqtag = line1[8].split("_")[0]
        mutations.append(refaa + start + line1[9])
        mut_nb += 1
    diffseq_annot = [conseqtag, mutations, starts, ends, mut_nb]
    return(diffseq_annot)

def get_barcode_cluster(cdhit_bak_cluster_file):
    with open(cdhit_bak_cluster_file) as f:
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

def get_mut_annot(mfa_to_annotate, reference_fa_file, barcode_cluster_dic):
    fasta_sequences = SeqIO.parse(open(mfa_to_annotate),'fasta')
    nb_fa_records = len(list(SeqIO.parse(open(mfa_to_annotate),'fasta')))
    ref_fasta_seq = str(SeqIO.read(open(reference_fa_file),'fasta').seq).upper()
    counter = 0
    cluster_number = 0
    dict_annot = {}
    with tempfile.NamedTemporaryFile(mode = "w") as my_fasta:
        for fasta in fasta_sequences:
            fabarcode = fasta.id
            fabarcode_clean = fabarcode.split("_")[0] # to remove _{d}
            counter += 1
            if ref_fasta_seq.upper() in str(fasta.seq).upper() or str(fasta.seq).upper() in ref_fasta_seq.upper(): # check if conseq is mutated if not annotate as WT and continue
                key = int(barcode_cluster_dic.get(fabarcode_clean)) # assign each annotation to unique cluster nb from cd-hit
                dict_annot[key] = [fabarcode_clean, "WT", "0", "0", "0"]
                continue
            SeqIO.write(fasta, my_fasta, "fasta")
            my_fasta.flush()
            output = subprocess.check_output(['diffseq', '-asequence'] +
                                             [str(reference_fa_file)] +
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
            key = int(barcode_cluster_dic.get(fabarcode_clean)) # assign each annotation to unique cluster nb from cd-hit
            dict_annot[key] = parse_diffseq_output(output, reference_fa_file)
            cluster_number += 1
    return dict_annot

def translate_mfasta(mfasta_file, outfile, frame):
    translated_records = []
    for sr in SeqIO.parse(mfasta_file, "fasta"):
        if frame == "1":
            sr.seq = sr.seq.translate(table="Bacterial")
            translated_records.append(sr)
        elif frame == "2":
            sr.seq = sr.seq[1:].translate(table="Bacterial")
            translated_records.append(sr)
        elif frame == "3":
            sr.seq = sr.seq[2:].translate(table="Bacterial")
            translated_records.append(sr)
        else:
            sys.exit("Please set reading frame 1, 2 or 3 or remove --readingframe option to use getorf for translation")
    SeqIO.write(translated_records, outfile, "fasta")

def main():
    
    # create outputdir if it doesn't exit
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    # write reference files
    if not os.path.exists(outfolder + "/reference"):
        os.makedirs(outfolder + "/reference")
    reference_nuc_file = outfolder + "/reference/nuc.fa"
    reference_prot_file = outfolder + "/reference/prot.fa"
    ref_nuc = SeqIO.read(os.path.abspath(args.ref), "fasta")
    with open(reference_nuc_file, "w") as output_handle:
    	SeqIO.write(ref_nuc, output_handle, "fasta")
    ref_prot = ref_nuc
    ref_prot.seq = ref_prot.seq.translate(table="Bacterial")
    with open(reference_prot_file, "w") as output_handle:
    	SeqIO.write(ref_prot, output_handle, "fasta")

    # run cd-hit-est to cluster identical amplicons
    print("\n" + "oooooooo Finding identical amplicons")
    logfile = open(outfolder + "/amplicons.log","a")
    cdhit = subprocess.call(["cd-hit-est", "-bak", "1", "-c", "1", "-G", "0", "-aL", "1", "-i", outfolder + "/amplicons.fna", "-o", outfolder + "/amplicons.cdhit.fna"], stdout=logfile)

    # run EMBOSS getORF to find ORF of each conseq or use Biopython translation if a reading frame is set
    if readingframe == "getorf":
        get_ORF_cdhit = getorf["-sequence", outfolder + "/amplicons.cdhit.fna", "-outseq", outfolder + "/amplicons.cdhit.faa", "-table", "11", "-minsize", ORFminsize, "-reverse", "FALSE"]
        get_ORF_cdhit()
        get_ORF_all = getorf["-sequence", outfolder + "/amplicons.fna", "-outseq", outfolder + "/amplicons.faa", "-table", "11", "-minsize", ORFminsize, "-reverse", "FALSE"]
        get_ORF_all()
    else:
        translate_mfasta(outfolder + "/amplicons.cdhit.fna",
                         outfolder + "/amplicons.cdhit.faa", readingframe)
        translate_mfasta(outfolder + "/amplicons.fna",
                         outfolder + "/amplicons.faa", readingframe)
        
    # get dictionnary of barcode for cd-hit clusters and all amplicons (protein and nucleotide sequences)
    barcode_cluster = get_barcode_cluster(outfolder + "/amplicons.cdhit.fna.bak.clstr")
    barcode_prot = get_barcode_fasta(outfolder + "/amplicons.cdhit.faa")
    barcode_nuc = get_barcode_fasta(outfolder + "/amplicons.cdhit.fna")
    barcode_all_prot = get_barcode_fasta(outfolder + "/amplicons.faa")
    barcode_all_nuc = get_barcode_fasta(outfolder + "/amplicons.fna")

    # run EMBOSS diffseq on each conseq.orf entry
    try:
        rm[outfolder + "/amplicons.effect"]()
    except:
        pass
    
    print("\n" +"oooooooo Annotating mutations")
    clust_annot_prot = get_mut_annot(outfolder + "/amplicons.cdhit.faa", reference_prot_file, barcode_cluster)
    clust_annot_nuc = get_mut_annot(outfolder + "/amplicons.cdhit.fna", reference_nuc_file, barcode_cluster)
    
    # write annotation files
    print("\n" + "oooooooo Writing annotation files")
    out_table = open(outfolder+"/amplicons.effect", "a") 
    out_table.write("\t".join(("barcode",
                              "sample",
                              "prot_mutation",
                              "prot_start",
                              "prot_end",
                              "nuc_mutation",
                              "nuc_start",
                              "nuc_end",
                              "prot",
                              "dna",
                               "reference_barcode",
                               "\n")))
    
    for barcode in barcode_all_nuc:
            
        nuc_annot = clust_annot_nuc.get(int(barcode_cluster.get(str(barcode))))
        if "WT" in nuc_annot:
            nuc_mut = "WT"
        else:
            nuc_mut = ",".join(nuc_annot[1])
        nuc_start = ",".join(nuc_annot[2])
        nuc_end = ",".join(nuc_annot[3])
        nuc = barcode_all_nuc.get(str(barcode))

        prot_annot = clust_annot_prot.get(int(barcode_cluster.get(str(barcode))))
        if prot_annot == None:
            prot_mut = "NA"
            prot_start = "NA"
            prot_end = "NA"
            prot = "NA"
        elif nuc_mut == "WT":
            prot_mut = "WT"
            prot_start = "0"
            prot_end = "0"
            prot = "NA"
        elif "NA" in prot_annot:
            prot_mut = "NA"
            prot_start = "NA"                                                                                                                                                                                                         
            prot_end = "NA"                                                                                                                                                                                                           
            prot = "NA"
        else:
            prot_mut = ",".join(prot_annot[1])            
            prot_start = ",".join(prot_annot[2])
            prot_end = ",".join(prot_annot[3])
            prot = barcode_all_prot.get(str(barcode))
                    
        ref_barcode = nuc_annot[0]

        annotation = "\t".join((barcode, sample_name, prot_mut, prot_start, prot_end, nuc_mut, nuc_start, nuc_end, prot, nuc, ref_barcode, "\n"))
        out_table.write(annotation)    

            
    out_table.close()
    os.remove("seq1out")
    os.remove("seq2out")

if __name__ == '__main__':
    main()
