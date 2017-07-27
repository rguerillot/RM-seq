# AMR-seq

# Usage

```
AMR-seq.pl [options] --R1 R1.fq.gz --R2 R2.fq.gz --refnuc FASTA --refprot FASTA --outdir DIR --barlen NN
  --help          This help.
  --debug!        Debug info (default '0').
  --R1=s          Read 1 FASTQ (default '').
  --R2=s          Read 2 FASTQ (default '').
  --refnuc=s      Reference gene that will be used for premapping filtering (fasta) (default '').
  --refprot=s     Reference protein that will be use for annotating variants (fasta) (default '').
  --outdir=s      Output folder (default '').
  --force!        Force overwite of existing (default '0').
  --barlen=i      Length of barcode (default '16').
  --minfreq=i     Minimum barcode frequency to keep (default '5').
  --cpus=i        Number of CPUs to use (default '1').
  --minsize=i     Minimum ORF size in bp used when annotating variants (default '200').
  --wsize=i       Word-size option to pass to diffseq for comparison with reference sequence (default '5').
  --subsample=i   Only examine this many reads (default '0').
  --keepfiles!    Do not delete intermediate files (default '0').
```
