#!/usr/bin/perl
use strict;
use warnings;
# coded by Qing Li on 07/29/16.
#Copyright (c) 2016 yandell lab. All rights reserved.


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "cleanUpReads.pl <leftReads.fq> <rightReads.fq> <adaptor.fa> outName\nneed to have bbmap, fqtrim,prinseq-lite in your path, the input reads must be named as xx.fq.gz or xx.txt.gz or xx.fastq.gz\n";

my$leftReads=$ARGV[0];
my$rightReads=$ARGV[1];
my$adaptor=$ARGV[2];
my$sample=$ARGV[3];
die $usage unless $sample;
die "wrong input file" unless $leftReads=~ /\.(fq|txt|fastq)\.gz/ && $rightReads=~ /\.(fq|txt|fastq)\.gz/;
#kmer frequency histogram and associated stats of raw reads:
system("reformat.sh in1=$leftReads in2=$rightReads out=$sample.interleaved.fastq.gz");
system("khist.sh in=$sample.interleaved.fastq.gz hist=khist.$sample.interleaved.txt -Xmx20g");
#adapter trimming of raw reads with fqtrim, may need generate a unique adapter file for each sample with the specific index sequences incorporated into the adapter, :
system("fqtrim -o fqtrim.fastq -f $adaptor -A -R -p 30 -l 80 $leftReads,$rightReads");
#quality trimming and filtering of reads with prinseq:
my($sample1)=$leftReads=~ /(\S+)\.(fq|txt|fastq)\.gz/;
my($sample2)=$rightReads=~ /(\S+)\.(fq|txt|fastq)\.gz/;
#this is to see the progress ...
system("gunzip $leftReads");
system("gunzip $rightReads");
system("prinseq-lite.pl -fastq $sample1.fqtrim.fastq -fastq2 $sample2.fqtrim.fastq -trim_left 1 -trim_qual_right 20 -trim_qual_type min -trim_qual_window 8 -trim_qual_step 1 -ns_max_n 0 -min_qual_score 10 -min_len 80 -out_good $sample.fqtrim.prinseq -out_bad null");

#error correction with ecc.sh and generation of kmer frequency histogram and associated stats of cleaned reads:
my$in1="$sample.fqtrim.prinseq"."_1.fastq";
my$in2="$sample.fqtrim.prinseq"."_2.fastq";
system("reformat.sh in1=$in1 in2=$in2 out=$sample.fqtrim.prinseq.interleaved.fastq.gz");
system("ecc.sh in=$sample.fqtrim.prinseq.interleaved.fastq.gz out=$sample.fqtrim.prinseq.ecc.interleaved.fastq.gz  histout=khist.$sample.fqtrim.prinseq.ecc.interleaved.txt -Xmx20g");
system("reformat.sh in=$sample.fqtrim.prinseq.ecc.interleaved.fastq.gz out1=$sample.fqtrim.prinseq.ecc.1.fastq out2=$sample.fqtrim.prinseq.ecc.2.fastq");
#---------------------------------- SUBS -------------------------------------
