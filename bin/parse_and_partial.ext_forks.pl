#!/usr/bin/perl
use forks;
use forks::shared;
use strict;
use warnings;
use Bio::DB::Fasta;


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "parse_and_partial.ext.pl class.out cpus mismatch_allowed\n";
#1 parse classifier output to print read id per transcript id;
#2 run transcriptome data mining script to select, deconcatenate, trinity assembly, blastx, annotate;
#3 since trinity run can only accomodate 1-2/run, when bash run this, only run 1-2 sp/run, so bx, annotate should give high cpus/count;
#4modify sub_contig id
#5pair the subcontig with original contig to run clustalw and cp the *cl.txt with limited # of mismatches to a new folder
my $class_out=$ARGV[0];
my $cpus=$ARGV[1];
my $misma_num=$ARGV[2];
die $usage unless $cpus;
my ($sp)=$class_out=~ /(\S+)\.tentative/;
system("mkdir $sp");
my%id;
open FH, $class_out;
while (<FH>) {
    next unless $_=~ /^C/;
    my@class=split /\t/, $_;
    if ($class[4]==2) {
	push @{$id{$class[3]}}, $class[1];
    }
}
close FH;
foreach my$rid(keys %id){
	open OUT, ">$sp/$rid";
	print OUT join("\n",@{$id{$rid}});
}
system("ls -1 $sp/ > $sp.id");
my @cid;
open(ID, "$sp.id");
while(my $id = <ID>){
                chomp($id);
		push @cid,$id;
}
my ($conca_reads)=</data1/qli/conus_transcriptome/concatenated/$sp.fa>;
my ($bx_db)=</archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein>;

#launch forks
#map here is to return a list that each element of the original list $_ becomes a list ref of 3 elements [$_, $sp, $conca_reads- the area between the {} is the code/change result block.
run_in_forks(\&select_reads, $cpus, map {[$_, $sp, $conca_reads]} @cid);
run_in_forks(\&deconcatenate_reads, $cpus, map {[$_, $sp]} @cid);
run_in_forks(\&trinity_assembly, ($cpus > 2) ? 2 : $cpus, map {[$_, $sp]} @cid);
run_in_forks(\&blastx, $cpus, map {[$_, $sp, $bx_db]} @cid);
run_in_forks(\&annotation, $cpus, map {[$_, $sp]} @cid);

#system("generate.bash_run_for_transcriptome_data_mining.pl --select_reads $sp.id $sp /data2/qli/tax_run/conus_trimmed_RNAseq_reads/concatenated/$sp.fa /archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein > $sp.ft");
#system("bash -x $sp.ft");#need be in forks
#system("generate.bash_run_for_transcriptome_data_mining.pl --deconcatenate_reads $sp.id $sp /data2/qli/tax_run/conus_trimmed_RNAseq_reads/concatenated/$sp.fa /archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein > $sp.dc");
#system("bash -x $sp.dc");#need be in forks
#system("generate.bash_run_for_transcriptome_data_mining.pl --trinity_assembly $sp.id $sp /data2/qli/tax_run/conus_trimmed_RNAseq_reads/concatenated/$sp.fa /archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein > $sp.tri");
#system("bash -x $sp.tri");#need be in forks

#system("generate.bash_run_for_transcriptome_data_mining.pl --blastx $sp.id $sp /data2/qli/tax_run/conus_trimmed_RNAseq_reads/concatenated/$sp.fa /archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein > $sp.bx");
#system("bash -x $sp.bx");#need be in forks
#system("generate.bash_run_for_transcriptome_data_mining.pl --annotation $sp.id $sp /data2/qli/tax_run/conus_trimmed_RNAseq_reads/concatenated/$sp.fa /archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein > $sp.anot");
#system("bash -x $sp.anot");#need be in forks
##modify subcontig id
open OUT1, ">$sp.partial.ext.fa";
my@fa_name_files=<./$sp/trin/annotated.*fa.txt>;
foreach my $file(@fa_name_files){
	my($file_id)=$file=~ /annotated.(\S+).fa.txt/;
	open FH, $file;
	while(<FH>) {
                        if (/^>/) {
                                $_=~ s/>/>$file_id._s/;
                                print OUT1;
                        }else{
                                print OUT1;
                        }
        }
        close FH;
}

##pair the original and partial ext to prep for clustalw
system("blastx -db /archive02/qli/data/cnserver_db/bx_db/ncbi_cn_unipro_sam_protein -query $sp.tentativeseq.contig.truncate.rf.fa -num_threads 20 -outfmt '6 std qframe sframe' -evalue 1e-3 -word_size 3 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -comp_based_stats t -seg no -soft_masking true -out $sp.ori.bx.out");
system("annotate.ncbi.pl $sp.ori.bx.out $sp.tentativeseq.contig.truncate.rf.fa > $sp.annotated.ori.fa.txt");
system("mkdir $sp/clustalw");
#system("prep_partial_ext_fa4clustalw.pl $sp.partial.ext.fa $sp.annotated.ori.fa.txt");
my ($ext_fa_file) = <$sp.partial.ext.fa>;
my ($ori_fa_file) = <$sp.annotated.ori.fa.txt>;
#my$ori_nt_fa_file=get_nt($ori_fa_file);
#my$ext_nt_fa_file=get_nt($ext_fa_file);
get_nt($ori_fa_file);
get_nt($ext_fa_file);
my ($ori_nt_fa_file)= <$ori_fa_file.nt>;
my ($ext_nt_fa_file)= <$ext_fa_file.nt>;

my($ori_id,$id,@anot);
my $count=0;
my $db=Bio::DB::Fasta->new($ori_nt_fa_file);
open (FH, $ext_nt_fa_file) or die "Can't open $ext_nt_fa_file for reading: $!\n";
while (defined (my$line=<FH>)){
    my $OUT;#define the scalar variable out of the scope
        if($line=~ />/){
                chomp ($line);
                ($id)=$line=~ />(\S+)/;
		($ori_id)=$line=~ />(\S+)\._s/;
                my$header_o=$db->header($ori_id);
                $line=~ s/>//;
                my$header_e=$line;
                my$anot_o=get_anot($header_o);
                my$anot_e=get_anot($header_e);
                push @anot,[$anot_o,$anot_e];
		my $seq=$db->get_Seq_by_id($ori_id);
                my $seqstr=$seq->seq;
                open $OUT,">$sp/clustalw/$count" or die "ERROR: Could not open filehandle to $count: $!";
                print $OUT ">$ori_id\n$seqstr\n>$id\n";
        }else{
        	open $OUT,">>$sp/clustalw/$count" or die "ERROR: Could not open filehandle to $count: $!";
        	print $OUT $line;
        	$count++;
        	close $OUT;
        }
}
close(FH);
open OUT2,">$sp/clustalw/$sp.annotation.table.txt";
system("mkdir $sp/clustalw/$sp.for_check/");
#run clustalw in forks
my @list;
my@indx=(0..$#anot);
foreach my $index (@indx){
    print OUT2 "$index\t";
    print OUT2 join("\t",@{$anot[$index]});
    print OUT2 "\n";

    if ($anot[$index]->[0] eq $anot[$index]->[1]){
	system("cp $sp/clustalw/$index $sp/clustalw/$sp.for_check/");
	push(@list, [$index]);
    }
}
run_in_forks(\&clustalw, $cpus, @list);

system("mkdir $sp/clustalw/$sp.for_check/$sp.$misma_num.mismatch");
my @cl_alignment_name_files=<$sp/clustalw/$sp.for_check/*.cl.txt>;
foreach my $cl_file(@cl_alignment_name_files){
	open FH3, $cl_file;
	my %align1;
        my %align2;
        while(defined(my $line1=<FH3>)){
                chomp($line1);
                my@seq1=split /\s+/,$line1;
                if(@seq1==2 && $line1 =~ /\w+/){
                        my $line2=<FH3>;
                        chomp($line2);
                        my@seq2=split /\s+/,$line2;
                        $align1{$seq1[0]}.=$seq1[1];
                        $align2{$seq2[0]}.=$seq2[1];
		}
	}
	my ($id1)=keys %align1;
        my ($id2)=keys %align2;
        my @array1=split //,$align1{$id1};
        my @array2=split //,$align2{$id2};
        my $gap=0;
        my $misma=0;
        for (my $j=0;$j<@array1;$j++){
		if($array2[$j] eq "-" || $array1[$j] eq "-"){
                	$gap++;
                	next;
		}elsif($array1[$j] ne $array2[$j]){
                	$misma++;
		}
	}
	if ($misma<= $misma_num){
		my ($fa_file_id)=$cl_file=~ /for_check\/(\d+)\.cl.txt/;
		system("cp $sp/clustalw/$sp.for_check/$fa_file_id $sp/clustalw/$sp.for_check/$sp.$misma_num.mismatch/");
		system("cp $cl_file $sp/clustalw/$sp.for_check/$sp.$misma_num.mismatch/");
	}
}

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
sub get_anot {
        my $header=shift;
        my $annotation;
        my @annot=split(/\t/,$header);
        my @note=split(/:/,$annot[1]);
        if ($note[1]=~ /OS=/){
                my @swp=split(/OS=/,$note[1]);
                $annotation=$swp[0];
        }else{
                $annotation=$note[1];
        }
        return $annotation;
}



sub get_nt {
	my$anot_fa_file=shift;
 	open(FH1,$anot_fa_file);
	open(OUT1,">$anot_fa_file.nt");
	while(defined(my$line1=<FH1>)){
        	if($line1=~ />/){
                	print OUT1 $line1;	
			my$line2=<FH1>;
			$line2=~ s/^nucleotide://;
			print OUT1 $line2;
		}
	}
#	my ($anot_nt_fa) = <$anot_fa_file.nt>;
#	return $anot_nt_fa;
 }

sub run_in_forks{	
    my $code_ref = shift;
    my $cpus = shift;
    
    #the rest of the @_ are @args 
    while(@_){
	if($cpus == 1){
	    my $args = shift @_;
	    &$code_ref(@$args);# if no threads, then directly dereference the code reference and run it
	}
	elsif(threads->list(threads::running) < $cpus){
	    my $args = shift @_;
	    threads->new($code_ref, @$args);#if some are finished then add new threads
	}
	elsif(my ($thr) = threads->list(threads::joinable)){
	    $thr->join;#join those finished not joined yet threads
	}
	else{
	    sleep 0.1;
	}
    }
    $_->join foreach(threads->list);#block until all threads finish running
}

sub select_reads{
	my $id=shift;
	my $id_dir=shift;
	my $conca_reads=shift;
	system("xargs samtools faidx $conca_reads < $id_dir/$id > $id_dir/$id.fa");
}
sub deconcatenate_reads{
	my $id=shift;
        my $id_dir=shift;
	system("~/anaconda/bin/python ~/tool/Taxonomer_py/utilities/deconcatenate_fasta.py -i  $id_dir/$id.fa -o $id_dir/$id");
}

sub trinity_assembly{
	my $id=shift;
        my $id_dir=shift;
	system("Trinity --seqType fa --max_memory 30G --left $id_dir/$id"."_1.fa --right $id_dir/$id"."_2.fa --CPU 6 --bflyHeapSpaceMax 5G --bflyCPU 2 --KMER_SIZE 32 --min_contig_length 75 --min_per_id_same_path 99 --max_diffs_same_path 1 --max_internal_gap_same_path 3 --output $id_dir/trin/$id.trinity --full_cleanup");
}
sub blastx{
	my $id=shift;
        my $id_dir=shift;
	my $bx_db=shift;
	system("blastx -db $bx_db -query $id_dir/trin/$id.trinity.Trinity.fasta -num_threads 1 -outfmt '6 std qframe sframe' -evalue 1e-3 -word_size 3 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -comp_based_stats t -seg no -soft_masking true -out $id_dir/trin/$id.bx.out");
}
sub annotation{
	 my $id=shift;
         my $id_dir=shift;
	 system("annotate.ncbi.pl $id_dir/trin/$id.bx.out $id_dir/trin/$id.trinity.Trinity.fasta > $id_dir/trin/annotated.$id.fa.txt");
}

sub clustalw{
    my $i = shift;
    system("clustalw2 -INFILE=$sp/clustalw/$sp.for_check/$i -align -type=DNA -outfile=$sp/clustalw/$sp.for_check/$i.cl.txt -output=clustalw");
    return;             
}

