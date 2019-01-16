#!/usr/bin/env perl

use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);

sub natatime ($@)
{
    my $n = shift;
    my @list = @_;

    return sub
    {
        return splice @list, 0, $n;
    }
}


my(@Options, $debug, $R1, $R2, $barlen, $outdir, $basequal, $refnuc, $minfreq, $force, $cpus, $minsize, $wsize, $subsample, $keepfiles, $noaln);
setOptions();

# Options
-r $R1 or err("Can't read --R1 $R1");
-r $R2 or err("Can't read --R2 $R2");
$barlen > 0 or err("Please provide length of barcode with --barlen");
$outdir or err("Please provide output folder with --outdir");

# Make output folder
if (-d $outdir) {
  if ($force) {
    msg("Forced removal of existing --outdir $outdir (please wait)");
    remove_tree($outdir);
  }
  else {
    err("Folder '$outdir' already exists. Try using --force");
  }
}
make_path($outdir);

# Subsampling reads if --subsample is set 
if ($subsample > 0) {
	msg("\noooooooo Subsampling $subsample reads");
	run_cmd("seqtk sample -s100 $R1 $subsample > $outdir/$subsample-reads-R1.fq");
	run_cmd("seqtk sample -s100 $R2 $subsample > $outdir/$subsample-reads-R2.fq");
	$R1 = ("$outdir/$subsample-reads-R1.fq");
	$R2 = ("$outdir/$subsample-reads-R2.fq");
}

# Trim read from 3' end if base quality below threshold
msg("\noooooooo Trimming end of reads when base quality is below $basequal");
run_cmd("trimmomatic PE -threads $cpus -phred33 \Q$R1\E \Q$R2\E $outdir/trimmed-R1.fq $outdir/trimmed-R1-unpaired.fq $outdir/trimmed-R2.fq $outdir/trimmed-R2-unpaired.fq TRAILING:$basequal > /dev/null 2>&1");

# Keep reads that map on reference
msg("\noooooooo Keeping reads mapping to reference");
run_cmd("mkdir $outdir/reference");
run_cmd("cp $refnuc $outdir/reference/ref.fa");
run_cmd("bwa index $outdir/reference/ref.fa >> $outdir/amplicons.log 2>&1");
run_cmd("bwa mem -t $cpus $outdir/reference/ref.fa $outdir/trimmed-R1.fq $outdir/trimmed-R2.fq 1> $outdir/out.sam 2>> $outdir/amplicons.log");
run_cmd("samtools view -Sb $outdir/out.sam 1> $outdir/out.bam 2>> $outdir/amplicons.log");
run_cmd("samtools view -b -f 0x2 $outdir/out.bam 1> $outdir/mapped.bam 2>> $outdir/amplicons.log");
run_cmd("bamToFastq -i $outdir/mapped.bam -fq $outdir/mapped_reads_R1.fq -fq2 $outdir/mapped_reads_R2.fq >> $outdir/amplicons.log 2>&1");

# Overlap reads
msg("\noooooooo Overlapping reads");
run_cmd("pear -u 0 -v 20 -j $cpus -f $outdir/mapped_reads_R1.fq -r $outdir/mapped_reads_R2.fq -o $outdir/reads >> $outdir/amplicons.log 2>&1");

# Count
msg("\noooooooo Counting barcodes");
#my $head = $subsample > 0 ? " head -n $subsample | " : "";
run_cmd("cat \Q$outdir/reads.assembled.fastq\E | paste - - - - | cut -f 2 | cut -c1-$barlen | sort | uniq -c | sort -nr > \Q$outdir/amplicons.barcodes\E");

# Filter out low freq
my %keep;
my $found=0;
open COUNT, '<', "$outdir/amplicons.barcodes";
while (<COUNT>) {
  chomp;
  my($count, $barcode) = split ' ';  # special ' ' whitespacer
  $keep{$barcode}=$count if $count >= $minfreq;
#  msg("# $count $barcode");
  $found++;
}
my $kept = scalar keys %keep;
msg("Found $found barcodes, keeping $kept with frequency >= $minfreq");

# Go back and bin reads
msg("\noooooooo Binning reads and creating FASTA files");
my %seq;
open RAW, "-|", "cat \Q$outdir/reads.assembled.fastq\E | paste - - - - | cut -f 2";
while (my $dna = <RAW>) {
  chomp $dna;
  my $barcode = substr $dna, 0, $barlen;
  next unless $keep{$barcode};
  push @{ $seq{$barcode} }, substr $dna, $barlen;
}

# Write FASTA files
my $counter=0;
for my $barcode (keys %keep) {
  print STDERR "\rWriting $barcode ", ++$counter, "/$kept";
  open my $fh, '>', "$outdir/$barcode.mfa";
  my $seqs = $seq{$barcode};
  my $nseq = scalar(@$seqs);
  for my $i (1 .. $nseq) {
    print $fh ">$barcode.$i\n", $seqs->[$i-1], "\n";
  }
  close $fh;
}

# Create consensus
msg("\n\noooooooo Creating consensus sequences (please wait)");
my $nbjobs = $cpus;
my $clustalo_cpu = 1;
if ($noaln != 0) {
	#skip clustalo alignment and run conseq only
	run_cmd("find $outdir -maxdepth 1 -name '*.mfa' | nice parallel --bar --progress -j $nbjobs \'cons -sequence {} -filter -name {/.} -plurality 0.8\' >> \Q$outdir\E/amplicons.fna\E", 1);
} else {
	#run clustalo alignment + conseq	
	run_cmd("find $outdir -maxdepth 1 -name '*.mfa' | nice parallel --bar --progress -j $nbjobs \'clustalo -i {} --outfmt=fa --threads=$clustalo_cpu | cons -filter -name {/.} -plurality 0.8\' >> \Q$outdir\E/amplicons.fna\E", 1);
}

# Cleanup
unless ($keepfiles) {
  msg("\noooooooo Deleting intermediate files");
  unlink "$outdir/$_.mfa", "$outdir/out.sam", "$outdir/out.bam", "$outdir/mapped.bam", "$outdir/trimmed-R1.fq", "$outdir/trimmed-R2.fq", "$outdir/trimmed-R1-unpaired.fq" , "$outdir/trimmed-R2-unpaired.fq", "$outdir/mapped_reads_R1.fq", "$outdir/mapped_reads_R2.fq" for keys %keep;
  unlink <$outdir/reads.*>;
}

exit(0);

#----------------------------------------------------------------------

sub run_cmd {
  my($cmd, $quiet) = @_;
  msg("Running: $cmd") unless $quiet;
  system($cmd)==0 or err("Error $? running command");
}

#----------------------------------------------------------------------
sub msg {
  print STDERR "@_\n";
}
      
#----------------------------------------------------------------------
sub err {
  msg(@_);
  exit(1);
}

#----------------------------------------------------------------------
sub num_cpus {
  my($num)= qx(getconf _NPROCESSORS_ONLN); # POSIX
  chomp $num;
  return $num || 1;
}
   
#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"R1=s",  VAR=>\$R1, DEFAULT=>'', DESC=>"Read 1 FASTQ"},
    {OPT=>"R2=s",  VAR=>\$R2, DEFAULT=>'', DESC=>"Read 2 FASTQ"},
    {OPT=>"refnuc=s",  VAR=>\$refnuc, DEFAULT=>'', DESC=>"Reference gene that will be used for premapping filtering (fasta)"},
    {OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'', DESC=>"Output folder"},
    {OPT=>"force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force overwite of existing"},
    {OPT=>"barlen=i",  VAR=>\$barlen, DEFAULT=>16, DESC=>"Length of barcode"},
    {OPT=>"minfreq=i",  VAR=>\$minfreq, DEFAULT=>5, DESC=>"Minimum barcode frequency to keep"},
    {OPT=>"basequal=i",  VAR=>\$basequal, DEFAULT=>30, DESC=>"Minimum base quality threshold used for trimming the end of reads (trimmomatic TRAILING argument)"},
    {OPT=>"cpus=i",  VAR=>\$cpus, DEFAULT=>&num_cpus(), DESC=>"Number of CPUs to use"},
    {OPT=>"minsize=i",  VAR=>\$minsize, DEFAULT=>200, DESC=>"Minimum ORF size in bp used when annotating variants"},
    {OPT=>"wsize=i",  VAR=>\$wsize, DEFAULT=>5, DESC=>"Word-size option to pass to diffseq for comparison with reference sequence"},
    {OPT=>"subsample=i",  VAR=>\$subsample, DEFAULT=>0, DESC=>"Only examine this many reads"},
    {OPT=>"keepfiles!",  VAR=>\$keepfiles, DEFAULT=>0, DESC=>"Do not delete intermediate files"},
	{OPT=>"noaln!",  VAR=>\$noaln, DEFAULT=>0, DESC=>"Skip reads alignment when generating consensus (to use for indel quantification only)"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] --R1 R1.fq.gz --R2 R2.fq.gz --refnuc FASTA --outdir DIR --barlen NN\n"; 
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
