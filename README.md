# RM-seq: Resistance Mutation SEQuencing

Analysis bioinformatic pipeline for high-throughput identification and quantification of large repertoires of resistance conferring mutations.
RM-seq is an amplicon-based, deep-sequencing technique founded on the single molecule barcoding technique. We have adapted this method in order to identify and quantify at high-throughput, the repertoire of mutations that confer resistance to a given antibiotic.
RM-seq allows to both correct sequenced read errors generated during sequencing and to accurately quantify mutations by correcting PCR amplification bias generated during sequencing library preparation. During the first step of amplicon library preparation, a linear PCR (primer extension) with a primer comprising a tail with degenerated bases (all possible bases) introduce a unique barcode to DNA template molecules. Therefore a barcode is assigned not just to all the molecules from a certain sample (indexing), but to all molecules being amplified and sequenced. RM-seq pipeline use these barcodes to generate an error-corrected consensus sequence of the initial template variant. Counting the barcodes indroduced before exponential amplification by PCR of the template allows to accurately quantify each genetic variants from genomic DNA extracted from complex population of resistant clones (eg. pools of 10,000 resistant colonies selected by an antibiotic in vitro).
A complete descrition of the RM-seq method will be available soon (article submitted)

## Is this the right tool for me?

1. To be able to us this pipeline you need to have sequenced amplicon library with molecular barcodes.
2. It only supports paired-end FASTQ reads (including .gz compressed fastq files).
3. It needs paired reads that are overlapping.
3. It needs bwa aligner and EMBOSS to be installed.
4. It needs a reference fasta sequences of the sequenced gene (DNA and protein sequence).
4. It's written in Python and Perl.

## Installation

### Install RM-seq pipeline
```
pip3 install rmseq
```
 
### Dependencies
RM-seq has the following package dependencies:
* EMBOSS >= 6.6 for `clustalo`, `cons`, `getorf`, `diffseq`
* bwa >= 0.7.15
* samtools >= 1.3
* pear >= 0.9.10
* cd-hit >= 4.7
* trimmomatic >= 0.36
* python modules: `plumbum`, `Biopython`

If you are using the [OSX Brew](http://brew.sh/) or [LinuxBrew](http://linuxbrew.sh/) packaging system:
```
brew tap homebrew/science
brew tap tseemann/bioinformatics-linux
brew install bwa
brew install samtools
brew install pear
brew install cd-hit
brew install trimmomatic
pip3 install plumbum
```

## Usage

### Quick start
```
rmseq run -1 <read1> -2 <read2> -n <ref_fasta_nuc> -p <ref_fasta_prot> -o <outdir>
```

### Get help
```
rmseq run -h
    usage: rmseq run [options]

    Run the pipeline

    optional arguments:
      -h, --help            show this help message and exit
      -d, --debug_on        Switch on debug mode.
      -1 R1, --R1 R1        Path to read pair 1
      -2 R2, --R2 R2        Path to read pair 2
      -n REFNUC, --refnuc REFNUC
                            Reference gene that will be used for premapping
                            filtering (fasta).
      -p REFPROT, --refprot REFPROT
                            Reference protein that will be use for annotating
                            variants (fasta).
      -o OUTDIR, --outdir OUTDIR
                            Output directory.
      -f, --force           Force overwite of existing.
      -b BARLEN, --barlen BARLEN
                            Length of barcode (default 16)
      -m MINFREQ, --minfreq MINFREQ
                            Minimum barcode frequency to keep (default 5)
      -c CPUS, --cpus CPUS  Number of CPUs to use (default 72)
      -r MINSIZE, --minsize MINSIZE
                            Minimum ORF size in bp used when annotating variants
                            (default 200)
      -w WSIZE, --wsize WSIZE
                            Word-size option to pass to diffseq for comparison
                            with reference sequence (default 5)
      -s SUBSAMPLE, --subsample SUBSAMPLE
                            Only examine this many reads.
      -k, --keepfiles       Keep the intermediate files. Default is to remove
                            intermediate files

## Output

RM-seq produces a tap-separated output file called amplicons.effect with the following columns:

Column | Example | Description
-------|---------|------------
barcode | GACACAACTGAGATTA | The sequence of the barcode
sample | Rifampicin1 | The output folder name
aa_mutation | H481N | The annotation of the amino acid change (Histidine residue 481 substituted by Asparagine)
start |  481 | start coordinate of the mutation 
end | 481 | end coordinate of the mutation
orf | VRPPDKNNRFVGLYCTLVSCVEWLAAAVVLYFCGVIVDAHVSFMSFIAIFIIAALSGLVSFIPGGFGAFDLVVLLGFKTLGVPEEKVLLMLLLYRF | the protein sequence of the consensus sequence
dna | GGTTAGACCACCCGATAAAAACAATCGTTTTGTAGGATTGTACTGCACTTTAGTGTCGTGTGTTGAATGGTTAGCAGCTGCAGTTGTATTATATTTCTGTGGTGTAATTGTTGACGCTCATGTATCATTCATGTCCTTTATTGCAATATTTATCATTGCTGCATTATCAGGTTTAGTCAGCTTTATTCCTGGTGGTTTCGGCGCTTTCGATTTAGTTGTATTACTAGGATTTAAAACTTTAGGTGTCCCTGAGGAAAAAGTATTATTAAT | The dna sequence of the consensus sequence

The other files produced by RM-seq are:

File name | Description
----------|------------
amplicons.nuc | Multifasta file containing all the consensus nucleotide sequence (header of sequence is the barcode)
amplicons.orf | Multifasta file containing all the consensus protein sequence (header of sequence is the barcode)
amplicons.barcodes | Table with the count of each barcode sequence
amplicons.cdhit | Multifasta file containing all the unique consensus nucleotide sequence (header of sequence is the barcode)


## Issues

Please report problems to the [Issues Page](https://github.com/rguerillot/RM-seq/issues).

## Author

Romain Guerillot | Torsten Seemann | Mark Schultz

