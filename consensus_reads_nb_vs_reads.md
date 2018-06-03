## Subsample reads file and run rmseq

### Set variables in your bash terminal

```
sampled_reads=0
increment=50000
total_reads_number=1000000 # do not exceed the number of reads in your fastq file
fastq_file_R1=my_R1.fq # replace with path to your R1 fastq file
fastq_file_R2=my_R2.fq # replace with path to your R2 fastq file
refnuc=my_ref_nuc.fa # replace with path to your reference nucleotide fasta file
refprot=ref_prot.fa # replace with path to your reference protein fasta file
```

### Subsample reads with seqtk and run rmseq

```
while [ "$sampled_reads" -lt "$total_reads_number" ]
do
  sampled_reads=$((sampled_reads+increment))
  echo "oooo Subsampling $sampled_reads reads with seqtk oooo"
  seqtk sample -s100 $fastq_file_R1 $sampled_reads > $sampled_reads_R1.fq
  seqtk sample -s100 $fastq_file_R2 $sampled_reads > $sampled_reads_R2.fq
  echo "oooo Running rmseq on $sampled_reads subsampled reads"
  rmseq run $sampled_reads_R1.fq $sampled_reads_R2.fq $refnuc $refprot $sampled_reads
done
```

### Concatenate all amplicoons.effect file obtained from different size of subsampling 

```
cat ./*_rmseq_outdir/amplicons.effect >> all.amplicons.tab
```

## Plot the number of consensus reads (from 10 reads) versus the number of reads (depth of sequencing) with Rstudio

### Install the tidyverse package
```
install.packages("tidyverse")
```

### Load tidyverse package
```
library(tidyverse)
```

### import all your consensus amplicon table into R
```
df_reads_subsampling <- read.table(file = "all_amplicons.tab", header = T)
```

### count the number of consensus amplicon obtained at each sequencing depth (fastq subsampling)

```
df_barcode_count <- df_reads_subsampling %>%
  group_by(sample) %>%
  count(sample, sort = T) %>%
  mutate(barcode_nb = n, reads_nb = as.integer(sample))
```

### plot consensus reads versus number of reads

```
ggplot() +
stat_smooth(data = df_barcode_count, aes(x=reads_nb, y = barcode_nb)) +
xlab(label = "number of reads")+
ylab(label = "number of conseq\n(consensus of 10 reads)")
```
