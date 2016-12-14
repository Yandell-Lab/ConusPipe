#!/usr/bin/perl 


use strict;
use warnings;
use lib '/home/qli/lib';
use Bio::DB::Fasta;
use Fasta_reader;
use Data::Dumper;
use Math::Random qw(random_negative_binomial);


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "simulatePipe.pl ref_transcriptomefasta_file totalTranscriptsNumber runfile referenceBlastDbName readLength insertSize outName\n
This script is to 
1.duplicate transcripts in transcriptome in a negative binormial distribution to mimic the real expression profile of transcrtips;
2. simulate reads from duplicated transcriptome;
3. separate simulated reads;
4. run serials of assmebly and bench mark for each assembly\n ";
my $refFa = $ARGV[0];
my $totalTransNum = $ARGV[1];
my $runFile=$ARGV[2];
my $rBnDbName=$ARGV[3];
my $readLen=$ARGV[4];
my $insertSize=$ARGV[5];
my $outName=$ARGV[6];
unless ($outName && -r $refFa){
    print $usage;
    exit;
}
#1.duplicate transcripts in transcriptome in a negative binormial distribution to mimic the real expression profile of transcrtips;
my @array_non_zero;
my @array = random_negative_binomial($totalTransNum, 1, 0.016);
foreach my $num (@array){
        $num=1 if $num==0;
        push @array_non_zero, $num;
}

my $i=0;
my $db=Bio::DB::Fasta->new($refFa);
open OUT,">$outName.dup.ref.fa";
open (FH, $refFa) or die "Can't open $refFa for reading: $!\n";
while (<FH>){
	if($_=~ />/){
		my ($id)= />(\S+)/;
		my $r_num=$array_non_zero[$i];
		dup_print($id,$r_num,\*OUT);
	$i++;	
	}
}
close(FH);
close(OUT);
#simulate reads from duplicated transcriptome;
system("cat $outName.dup.ref.fa | simLibrary -r $readLen -s same -x 100 -i $insertSize -b 1 > $outName.lib.fa");
system("cat $outName.lib.fa |simNGS -g 1e-3 -I --adapter -p paired $runFile > $outName.simulated_reads.fq 2>err.txt");
# separate simulated reads;
open(OUT1, ">$outName.1.fq");
open(OUT2, ">$outName.2.fq");
open(FH1, "<$outName.simulated_reads.fq");
while (defined(my $line1=<FH1>)){
        my$line2=<FH1>;
        my$line3=<FH1>;
        my$line4=<FH1>;
        print OUT1 "$line1.$line2.$line3.$line4" if $line1=~ /\/1/;
        print OUT2 "$line1.$line2.$line3.$line4" if $line1=~ /\/2/;
}
close(FH1);
close(OUT1);
close(OUT2);
#run serials of assmebly and bench mark for each assembly--from previous knowledge we know the 2 paramters need to test are kmer length and kmer coverage
my($kmerCov,$kmerLen);
system("makeblastdb -in $refFa -out $rBnDbName -dbtype nucl");
for $kmerCov(1..20){
	for $kmerLen(10..32){
		system("rm -rf $outName.kc.$kmerCov.kl.$kmerLen.Trinity");
		system("rm -rf $outName.kc.$kmerCov.kl.$kmerLen.Trinity.Trinity.fasta");
		system("Trinity --seqType fq --max_memory 30G --bypass_java_version_check --left $outName.1.fq --right $outName.2.fq --CPU 23 --bflyHeapSpaceMax 5G --bflyCPU 12 --KMER_SIZE $kmerLen --SS_lib_type RF --min_kmer_cov $$kmerCov --min_glue 10  --output $outName.kc.$kmerCov.kl.$kmerLen.trinity --full_cleanup");
system("makeblastdb -in $outName.kc.$kmerCov.kl.$kmerLen.Trinity.Trinity.fasta -out $outName.kc.$kmerCov.kl.$kmerLen -dbtype nucl");
		system("blastn -db $outName.kc.$kmerCov.kl.$kmerLen -query $refFa -num_threads 22 -outfmt 6 -evalue 1e-10 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out $outName.ref_vs_kc$kmerCov.kl$kmerLen.bn");
		system("blastn -db $rBnDbName -query $outName.kc.$kmerCov.kl.$kmerLen.Trinity.Trinity.fasta -num_threads 22 -outfmt 6 -evalue 1e-10 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out $outName.kc$kmerCov.kl$kmerLen.vs.refe.bn");
		my($o2tBn)=<$outName.ref_vs_kc$kmerCov.kl$kmerLen.bn>;	
		my($t2oBn)=<$outName.kc$kmerCov.kl$kmerLen.vs.refe.bn>;	
		my($ft)=<$outName.kc.$kmerCov.kl.$kmerLen.Trinity.Trinity.fasta>;
		reportSnSpTp($o2tBn,$t2oBn,$refFa,$ft,$outName);
	}

}


#---------------------------------------------subs---------------------------------------------------------
sub dup_print {
	my $id=shift;
	my $r_num=shift;
	my $fh=shift;
	my $header=$db->header($id);
	my $seq=$db->get_Seq_by_id($id);
	my $seqstr=$seq->seq;
	for(my $j=1;$j<=$r_num;$j++){
		print $fh ">$j"."_"."$header\n$seqstr\n";
	}
}

sub reportSnSpTp{
	my $o2t_blast_out =shift;
	my $t2o_blast_out =shift;
	my $fasta_file_o =shift;
	my $fasta_file_t =shift;
	my $output_prefix =shift;
	my $best_hit_o= best_hits($o2t_blast_out); # hash ref of blast output for original fasta as query .
	my $best_hit_t= best_hits($t2o_blast_out); # hash ref of blast output for trinity assembly fasta as query .
	my %rbh_o2t;#rbh pair original fasta id as key. 
	my %rbh_t2o;#rbh pair trinity assembly fasta id as key. 
	for my $o (keys %{$best_hit_o}){
		my $best_o=$best_hit_o->{$o}{db_id};
		my $best_t=$best_hit_t->{$best_o}{db_id};
		next unless ($o && $best_t);
		if ($o eq $best_t){
			$rbh_o2t{$o}=$best_o;
                $rbh_t2o{$best_o}=$o;
		}
	}
	## get sequence length info, coverage=query(ori seq) aligned length / query seq length, TP,FN,FP for sp and sn calculation
	my ($specif,$sensit,$ori_count);
           my $tp=0;
           my $fp=0;
           my $fn=0;
            my $fasta_reader = new Fasta_reader($fasta_file_o);
            open (my $ofh, ">$output_prefix.id_cov_rbh_qlength") or die $!;
            print $ofh join("\t", "#qseqid", "sseqid", "pident","pcover","qlen") . "\n";
            while (my $seq_obj = $fasta_reader->next()) {
                my $acc = $seq_obj->get_accession();
                $ori_count++;
                if (exists $rbh_o2t{$acc}) {
                    my $sequence = $seq_obj->get_sequence();
                    my $seq_length = length($sequence);
                    my $q_coverage= ($best_hit_o->{$acc}{query_match_len})/$seq_length;
                    print $ofh join("\t",$acc,$best_hit_o->{$acc}{db_id},$best_hit_o->{$acc}{percent_id},$q_coverage,$seq_length)."\n";
                    $tp=$tp+ $q_coverage;
                    my $test;
                }else{
                        $fn++;
                }
            }

            my $fasta_reader_B = new Fasta_reader($fasta_file_t);
	    while (my $seq_obj_B=$fasta_reader_B->next()){
                my $acc_B=$seq_obj_B->get_accession();
                $fp++ unless exists $rbh_t2o{$acc_B};
            }
	##  calculate sensitivity and specificity;
	$specif=$tp/($tp+$fp);
	$sensit=$tp/($tp+$fn);
	my$tp_rate=$tp/$ori_count;
	open ( FH,">$output_prefix.sn.sp.tpr.txt");
	print FH "$output_prefix\t$sensit\t$specif\t$tp_rate\n";

	close $ofh;
	close FH;




}			
sub best_hits{
	 my $file= shift;
    open (my $fh, $file) or die "Error, cannot open file $file";
    my %best_hits;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/, $_);
        my $query_id = $x[0];
        my $db_id = $x[1];
        my $percent_id = $x[2];
        my $query_start = $x[6];
        my $query_end = $x[7];
        my $db_start = $x[8];
        my $db_end = $x[9];

        my $evalue = $x[10];
        my $bitscore = $x[11];

        if ( (! exists $best_hits{$query_id}) || ($evalue < $best_hits{$query_id}->{evalue}) ) {

            $best_hits{$query_id} = {        db_id => $db_id,
                                             percent_id => $percent_id,
                                             query_start => $query_start,
                                             query_end => $query_end,
                                             db_start => $db_start,
                                             db_end => $db_end,
                                             evalue => $evalue,
                                             bitscore => $bitscore,

                                             query_match_len => abs($query_end - $query_start) + 1,
                                             db_match_len => abs($db_end - $db_start) + 1,

                                         };
        }
    }
    close $fh;
    return \%best_hits;
}
