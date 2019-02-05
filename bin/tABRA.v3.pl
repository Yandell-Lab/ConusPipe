#!/usr/bin/perl
use strict;
use warnings;
use Tie::IxHash;
use Bio::DB::Fasta;
# coded by Qing Li on 07/29/16.
#Copyright (c) 2016 yandell lab. All rights reserved.


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "tABRA.pl leftReads rightReads blastxDb pepLength bxCpus trinityCpus rsemCpus outName\n
#1  trinity assembly-optimized for long reads conotoxin assembly,label the ids with species, RSEM,E90 script(this annotation pipeline tABRA version2 matches trinity version 2.8.4 output matrix format) to generate statistic file(stderr to a output file), filter trinity assembly by tpm and isoformPercent(IsoPct) in RSEM.isoforms.results,get transcriptome statistics, blastx, annotate&get annotated ones' id;
#2 output intermediate file:label& sort anot.pep&trinity.filtered.fa with tpm.   
#3 remove anot ones from trinity assembly;
#4extract conotoxins from annotated ones, trim2longestORF(first Methionine/uppercase to the stop codon/end of seq, uppercase trim is prefered than methionine trim).remove redundant seq for trimmed ones and update tpm value;remove non-redundant seq without a stop codon, but still keep the non-redundant seq with stop codon even if they are without a methionine in the beginning.print out the anot.cono.nt.uniq.fa,anot.cono.pep.uniq.fa, initial statisic report is based on these 2 files. also print out the non-redundant seq with stop codon but without a methionine in the beginning as truncated seq for partial extention: anot.cono.nt.uniq.truncated.fa,anot.cono.pep.uniq.truncated.fa ;
#5.output the unannotated trinity seq,6frame translated pep seq in both chopped and whole length format (new!),       
## make a separate pipeline after this for potential novel conotixin discovery\n
#note: Sam requires stop codon to be in the end of the toxin seq;for nonCono, the stop codon was not removed either\ntABRA version2 matches trinity version 2.2.0 output matrix format\nif trinity run is already done, then need to put the output file in the run directory with name sample.trinity.Trinity.fasta, it will automatically skip trinity assembly\n";  

my$left_reads=$ARGV[0];
my$right_reads=$ARGV[1];
my$bx_db=$ARGV[2];
my$pep_RL=$ARGV[3];
my $bxCpus=$ARGV[4];
my $triCpus=$ARGV[5];
my $rsemCpus=$ARGV[6];
my $sample=$ARGV[7];
my$bflyCpus=$triCpus-5 if$triCpus;
die $usage unless $sample;
unless(-f "$sample.trinity.Trinity.fasta" ){
	system("Trinity --seqType fq --max_memory 30G --bypass_java_version_check --left $left_reads --right $right_reads --CPU $triCpus --bflyHeapSpaceMax 5G --bflyCPU $bflyCpus --KMER_SIZE 31 --SS_lib_type RF --min_kmer_cov 10 --min_glue 10  --output $sample.trinity --full_cleanup");
}
open FH,"$sample.trinity.Trinity.fasta" or die "$sample.trinity.Trinity.fasta does not exist!!";
open OUT,">$sample.trinity.Trinity.rf.fasta";
while(defined (my$line=<FH>)){
	$line=~ s/>/>$sample./;
	print OUT $line;
}
close FH;
close OUT;
system("align_and_estimate_abundance.pl --transcripts $sample.trinity.Trinity.rf.fasta --seqType fq --left $left_reads --right $right_reads --SS_lib_type RF --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir $sample.rsem_outdir --thread_count $rsemCpus ");
system("abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map $sample.trinity.Trinity.rf.fasta.gene_trans_map --out_prefix $sample.matrix.not_cross_norm $sample.rsem_outdir/RSEM.isoforms.results");
#note: if seq ids are modified, need to hack the script to get E90
system("contig_ExN50_statistic.pl $sample.matrix.not_cross_norm.isoform.TPM.not_cross_norm $sample.trinity.Trinity.rf.fasta | tee $sample.ExN50.stats");
#get transcriptome statistics
my($N50,$total,$assembNum)=calculateN50($sample);
my(%tpm,%isoPcL);
open FH,"$sample.rsem_outdir/RSEM.isoforms.results"; 
while(defined (my$line=<FH>)){
        next if $line=~ /^transcript_id/;
	chomp ($line);
        my@t=split /\t/,$line;
        $tpm{$t[0]}=$t[5];
        $isoPcL{$t[0]}=[$t[7],$t[2]];
}
close FH;
my($triFa)=<$sample.trinity.Trinity.rf.fasta>;
my$triDb=Bio::DB::Fasta->new($triFa);
open(OUT,">$sample.tri.nt.filtered.rsem.fa");
my@ids=sort{$tpm{$b} <=> $tpm{$a} } keys(%tpm);
foreach my $id (@ids) {
	my$tpm=$tpm{$id};
	my$isoPc=$isoPcL{$id}->[0];
	next if $tpm <=1 ||$isoPc <=1;
        my$header=">$id\tlen:$isoPcL{$id}->[1]\ttpm:$tpm";
        my$ntSeq=$triDb-> get_Seq_by_id($id);
        my $ntSeqstr  = $ntSeq-> seq;
        print OUT $header."\n$ntSeqstr\n";
}
close OUT;
system("blastx -db $bx_db -query $sample.tri.nt.filtered.rsem.fa -num_threads $bxCpus -outfmt '6 std qframe sframe' -evalue 1e-3 -word_size 3 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -comp_based_stats t -seg no -soft_masking true -max_target_seqs 1 -out $sample.bx.out");

my$nPHeaderSeq=annotation($sample);#{tid}->[entry,frame,ntSeq,pepSeq]
##for the annotated seq, separate into pep and nt, cono and nocono, get # of annotated seq
open OUT20, ">$sample.anot.transcriptome.pep.fa";
open OUT21, ">$sample.transcriptome.nt.fa";
open OUT1, ">$sample.anot.cono.nt.uniq.fa";
open OUT2, ">$sample.anot.cono.pep.uniq.fa";
open OUT8, ">$sample.anot.noCono.nt.uniq.fa";
open OUT9, ">$sample.anot.noCono.pep.uniq.fa";
my(%supFam,%pepTpm,$anotCount);
foreach my$tid(@ids){
        next unless defined $nPHeaderSeq->{$tid};
	if(@{$nPHeaderSeq->{$tid}}>1){
		$anotCount++;
		 my$imHeader=">$tid\tEntry:$nPHeaderSeq->{$tid}->[0]\tlen:$isoPcL{$tid}->[1]\ttpm:$tpm{$tid}\n";
		print OUT20 $imHeader.$nPHeaderSeq->{$tid}->[3]."\n";
		print OUT21 $imHeader.$nPHeaderSeq->{$tid}->[2]."\n";
		my$pepOrf=refinePep($nPHeaderSeq->{$tid}->[3]);
		if($nPHeaderSeq->{$tid}->[0]=~ /(conotoxin|Contulakin|contulakin|Contrypha|Conkunitzin|cokunitzin|contryphan|con-ikot|prepropeptide|insulin|prohormone-4|neuropeptide_L11|bombyxin|glucagon)/i && $nPHeaderSeq->{$tid}->[0]!~ /(receptor|binding|associated|induced|enzyme|enhancer|inhibitor)/i ){	
			if(not defined  $pepTpm{$pepOrf}){
                                $pepTpm{$pepOrf}=$tid ;
                        }elsif($isoPcL{$pepTpm{$pepOrf}}->[1]< $isoPcL{$tid}->[1]){
                                $tpm{$tid}+=$tpm{$pepTpm{$pepOrf}};
                                $pepTpm{$pepOrf}=$tid;
                        }else{
                                $tpm{$pepTpm{$pepOrf}}+=$tpm{$tid};                
                        }

		}else{
			#$pepOrf=~ s/\*//;
			my$ntOrf=trimNt($nPHeaderSeq->{$tid},$pepOrf);	
			print OUT8 $imHeader.$ntOrf."\n";				
			print OUT9 $imHeader.$pepOrf."\n";				

		}

	}else{          
                print OUT21 ">$tid\tEntry:noBxHit\tlen:$isoPcL{$tid}->[1]\ttpm:$tpm{$tid}\n$nPHeaderSeq->{$tid}->[0]\n";  
        }
}
close OUT20;
close OUT21;
close OUT8;
close OUT9;
open OUT5,">$sample.assembleStat.txt";
print OUT5 "N50:$N50\ntotal:$total\nassembly number:$assembNum\nbx hits:$anotCount\n"; 
close OUT5;
#new subroutine remove redundant and change tpm value(pass in %pepTpm & %tpm )
##Note:the value for %pepTpm should be the unmodified trinity transcript id!! 
my$uniqPT=removeRedundantUpdateTpm(\%pepTpm,\%tpm,\%isoPcL);
my%uniqPT;
foreach my$uPT(@$uniqPT){
                my$dbId=$pepTpm{$uPT};
                $uniqPT{$dbId}=$uPT;
}
open OUT3,">$sample.noTrim.anot.cono.nt.uniq.fa";
open OUT10, ">$sample.noTrim.anot.cono.truncated.nt.uniq.fa";
open OUT11, ">$sample.anot.cono.truncated.pep.uniq.fa";
@ids=sort{$tpm{$b} <=> $tpm{$a} } keys(%tpm);
foreach my $id (@ids) {
        next unless defined$uniqPT{$id};
        my($sFam)=$nPHeaderSeq->{$id}->[0]=~ /conotoxin\.(\S+)/;
	#print STDERR "$nPHeaderSeq->{$id}->[0]\t$sFam\n";#debug
	$supFam{$sFam}->[1]+=$tpm{$id};
	#print STDERR "$id\t$sFam\t$tpm{$id}\t$nPHeaderSeq->{$id}->[0]\t$supFam{$sFam}->[1]\n" if $sFam eq "B" ;#debug
	next unless $uniqPT{$id}=~ /\*/;
	#$uniqPT{$id}=~ s/\*//;
        $supFam{$sFam}->[0]++;
	 #print STDERR "$nPHeaderSeq->{$id}->[0]\t$supFam{$sFam}->[0]\n";#debug
	 my $mIndex=index($nPHeaderSeq->{$id}->[3],$uniqPT{$id});
        my$pepLen=length($uniqPT{$id});
        my $frame_strand = $nPHeaderSeq->{$id}->[1];
        my($frame)  = $frame_strand =~ /(\d+)/;
        my($strand) = $frame_strand =~ /([\+\-])/;
        $strand ||= '+';
        my $nt_start=$mIndex*3+$frame-1;
        my $ntLen=$pepLen*3;
        my $pep_orf=$uniqPT{$id};
        my$ntSeq=$nPHeaderSeq->{$id}->[2];
        my$ntFulLen=length($ntSeq);
        $ntSeq= revcom($ntSeq) if $strand eq '-';
        my $nt_orf  = substr $ntSeq, $nt_start, $ntLen;
         my$header="$id\tsupFam:$sFam\ttpm:$tpm{$id}";
         print OUT1 ">$header\tlen:$ntLen\n$nt_orf\n";
         print OUT2 ">$header\tlen:$pepLen\n$uniqPT{$id}\n";
         print OUT3 ">$header\tlen:$ntFulLen\n$ntSeq\n";
	if($pep_orf !~ /^[mM]/){
		print OUT11 ">$header\tlen:$pepLen\n$uniqPT{$id}\n";
		print OUT10 ">$header\tlen:$ntFulLen\n$ntSeq\n";
	}

}
close OUT1;
close OUT2;
close OUT3;
close OUT10;
close OUT11;

# get superfamily statistic( not including low tpm or low isoPct transcripts--filtered out,  but including truncated seq and complete seq with aligned part<10aa - not output )
open OUT4, ">$sample.bxAnot.supfam.stat.txt";
my@supFam=sort {$supFam{$b}->[1] <=> $supFam{$a}->[1]} keys(%supFam);
foreach my$sFam(@supFam){
        print OUT4 "$sFam\t";
        print OUT4 $supFam{$sFam}->[0] ||0;
        print OUT4 " \t$supFam{$sFam}->[1]\n";
}
close OUT4;


##remove the annotated seq from  total assembly.
my($tri_file)=<$sample.tri.nt.filtered.rsem.fa>;
open OUT6,">$sample.tri.no_anot.fa";
open OUT7,">$sample.tri.no_anot.pep.fa";
open OUT12,">$sample.tri.no_anot.whole.pep.fa";
my$db=Bio::DB::Fasta->new($tri_file);
open(FH, $tri_file) or die "Can't open $tri_file for reading: $!\n";
while (defined (my$line=<FH>)){
	if($line=~ /^>/){
		my($tid2)=$line=~ />(\S+)/;
		next if @{$nPHeaderSeq->{$tid2}}>1;
		my$seq=$db-> get_Seq_by_id($tid2);
		my $seqstr  = $seq-> seq;
		$line=~ s/>//;
		print OUT6 ">$line";
		print OUT6 "$seqstr\n";
		#run 6frame_translation in unannotated trinity assembly;
		for my $j (1..3){
                        translate_print_seq($tid2,$seqstr, $j,\*OUT7,$pep_RL);
                        translate_print_seq_whole($tid2,$seqstr, $j,\*OUT12);
                }
	}		

}
close OUT6;
close OUT7;
close OUT12;
close FH;

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------

sub trimNt{
	my$nPHS=shift;
	my$pepOrf=shift;
	my $mIndex=index($nPHS->[3],$pepOrf);
	my$pepLen=length($pepOrf);
	my $frame_strand = $nPHS->[1];
	my($frame)  = $frame_strand =~ /(\d+)/;
        my($strand) = $frame_strand =~ /([\+\-])/;
        $strand ||= '+';
        my $ntStart=$mIndex*3+$frame-1;
	my $ntLen=$pepLen*3;
	my$ntSeq=$nPHS->[2];
	my$ntOrf=substr $ntSeq, $ntStart, $ntLen;
	return $ntOrf;

}
sub refinePep{
	my$pep=shift;
	my @pepSeq;
	if($pep=~ /m[a-z]*[A-Z]+[A-Za-z]*\*/){
		@pepSeq=$pep=~ /(m[a-z]*[A-Z]+[A-Za-z]*\*)/g;
	}elsif($pep=~ /[A-Z]+[A-Za-z]*\*/){
		@pepSeq=$pep=~ /([A-Z]+[A-Za-z]*\*)/g;
	}elsif($pep=~ /m[a-z]*[A-Z]+[A-Za-z]*$/){
		@pepSeq=$pep=~ /(m[a-z]*[A-Z]+[A-Za-z]*)$/;
	}else{
		@pepSeq=$pep=~ /([A-Z]+[A-Za-z]*)$/;
	}
	# only keep fragment seq containing at least 5aa aligned part   
	my@upPepSeq;
	foreach(@pepSeq){
		push @upPepSeq,$_ if $_=~ /[A-Z]{5,}/;
	}
	next if scalar @upPepSeq==0; 
	my@sortedPepSeq=sort {length $b <=> length $a} @upPepSeq;	
	return $sortedPepSeq[0];
}

sub calculateN50{
	my$sample=shift;
	my($filename)=<$sample.trinity.Trinity.rf.fasta>;
	die "trinity file does not exist!" unless -r $filename;
	my $n;
my ($len, $total)=(0,0);
my @x;
open FH, $filename;
while (<FH>){
        if(/>/){
                if($len>0){
                        $total+=$len;
                        push @x,$len;
                 }
        $len=0;
        }else{
              s/\s//g;
              $len+=length($_);
        }
}
close FH;
#take care of the last seq
if ($len>0){
        $total+=$len;
        push @x,$len;
}
@x=sort{$b<=>$a} @x;
my ($count,$half,$N50)=(0,0);
for (my $j=0;$j<@x;$j++){
        $count+=$x[$j];
        if (($count>=$total/2)&&($half==0)){
                $n=@x;
                $N50=$x[$j];
                $half=$x[$j];
                last;
        }
}
return ($N50,$total,$n);
}
sub removeRedundant{
	my$pepSeq=shift;
	my@sorted=sort {length $b <=> length $a} @$pepSeq;
	my@uniqSeq;
	my$flag=0;
	foreach my $seq(@sorted){
		my$seqBL=uc($seq);
        	$flag=0;
        	foreach my$uSeq(@uniqSeq){
			my$uSeqBL=uc($uSeq);
                	if($uSeqBL=~ /$seqBL/){
                        	$flag=1;
                	}
        	}
        	if($flag==0){
                	push @uniqSeq,$seq;
        	}

	}
	return \@uniqSeq;
}



sub removeRedundantUpdateTpm{
        my$pepTpm=shift;
        my$tpm=shift;
        my$isoPcL=shift;
        my@pepTpm=keys %{$pepTpm};
        my@sorted=sort {length $b <=> length $a} @pepTpm;
        my@uniqSeq;
        my$flag=0;
        foreach my $seq(@sorted){
                my$pTid=$pepTpm->{$seq};
                my$seqBL=uc($seq);
                $seqBL=~ s/\*//g;
                $flag=0;
                foreach my$uSeq(@uniqSeq){
                        my$uTid=$pepTpm->{$uSeq};
                        my$uSeqBL=uc($uSeq);
                        $uSeqBL=~ s/\*//g;
                        if($uSeqBL eq $seqBL){
                                if($isoPcL{$uTid}->[1] >= $isoPcL{$pTid}->[1]){# here compared in original length
                                        $tpm->{$uTid}+=$tpm->{$pTid} if defined$tpm->{$pTid};
                                        delete $tpm->{$pTid};
                                        $flag=1;
                                }else{
                                        $tpm->{$pTid}+=$tpm->{$uTid} if defined$tpm->{$uTid};
                                        delete $tpm->{$uTid};
                                        my$index=0;
                                        $index++ until $uniqSeq[$index] eq $uSeq;
                                        splice(@uniqSeq,$index,1);
                                }
                        }elsif($uSeqBL=~ /$seqBL/){
                                $flag=1;
                                $tpm->{$uTid}+=$tpm->{$pTid} if defined$tpm->{$pTid};
                                 delete $tpm->{$pTid};
                        }
                }
                if($flag==0){
                        push @uniqSeq,$seq;
                }

        }
        return \@uniqSeq;
}



sub translate_print_seq {
        my $entry=shift;
        my $seq=shift;
        my $frame=shift;
	my $fh=shift;
	my $pep_RL=shift;
        my $pep = translate_with_frame($seq, $frame);
        my $rpep = translate_with_frame($seq, -$frame);
	print $fh $entry."_frame_$frame\n$pep\n";
	print $fh $entry."_frame_-$frame\n$rpep\n";
        my @peps=$pep=~ /([A-Z]{$pep_RL,})\*/g;
                for (my$i=0;$i<@peps;$i++){
                        print $fh ">$entry"."_frame_$frame"."_$i\n$peps[$i]\n";
                }

        my @rpeps=$rpep=~ /([A-Z]{$pep_RL,})\*/g;
                for (my$k=0;$k<@rpeps;$k++){
                        print $fh ">$entry"."_frame_-$frame"."_$k\n$rpeps[$k]\n";
                }
}
sub translate_print_seq_whole {
        my $entry=shift;
        my $seq=shift;
        my $frame=shift;
	my $fh=shift;
        my $pep = translate_with_frame($seq, $frame);
        my $rpep = translate_with_frame($seq, -$frame);
	print $fh ">$entry"."_frame_$frame\n$pep\n";
	print $fh ">$entry"."_frame_-$frame\n$rpep\n";


}


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
	close FH1;
	close OUT1;
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
	open FH, "$id_dir/$id.fa";
        open OUT1,">$id_dir/$id.1.fa";
        open OUT2,">$id_dir/$id.2.fa";
        open OUT3,">$id_dir/$id.s.fa";
	my($readFa)=<$id_dir/$id.fa>;
	my$db=Bio::DB::Fasta->new($readFa);
        my$rid;
	my$singletons=0;
        while(my$line=<FH>){
                chomp $line;
                if($line=~ />/){
                        ($rid)=$line=~ />(\S+)/;
			my$seq=$db-> get_Seq_by_id($rid);
                	my $seqstr  = $seq-> seq;
               		my$separator_indx=index($seqstr,'-');
			if($separator_indx==-1){
				print OUT3 ">$rid\n$seqstr\n";
				$singletons++;
			}else{	
				my$length=$separator_indx-1;	
				my$offset=$separator_indx+1;
 				my$read1=substr($seqstr,0,$length);
 				my$read2=substr($seqstr,$offset);
                        	print OUT1 ">$rid/1\n$read1\n";
                        	print OUT2 ">$rid/2\n$read2\n";
			}
                } 
        }
	close(FH);
	close(OUT1);
	close(OUT2);
	close(OUT3);
	system("rm $id_dir/$id.s.fa") if $singletons==0;
}

sub trinity_assembly{
	my $id=shift;
        my $id_dir=shift;
	system("Trinity --seqType fa --max_memory 30G --bypass_java_version_check --left $id_dir/$id".".1.fa --right $id_dir/$id".".2.fa --CPU 6 --bflyHeapSpaceMax 5G --bflyCPU 2 --KMER_SIZE 31 --SS_lib_type RF --min_kmer_cov 10 --min_glue 10 --min_contig_length 180 --output $id_dir/trin/$id.trinity --full_cleanup");
}
sub blastx{
	my $id=shift;
        my $id_dir=shift;
	my $bx_db=shift;
	system("blastx -db $bx_db -query $id_dir/trin/$id.trinity.Trinity.fasta -num_threads 1 -outfmt '6 std qframe sframe' -evalue 1e-4 -word_size 3 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -comp_based_stats t -seg no -soft_masking true -out $id_dir/trin/$id.bx.out");
}
sub annotation{
my $id=shift;
my $blast_file = "$id.bx.out";
my $t_file = "$id.tri.nt.filtered.rsem.fa";
my $e_thres = 1e-5;

open FH1,$blast_file;
my %toxin;
my %first;
while (<FH1>) {
        my @i = split /\t/, $_;

        if (not defined $toxin{$i[0]}{hit}) {
                $toxin{$i[0]}->{hit} = $i[1];
                $toxin{$i[0]}->{frame} = $i[12];
                $toxin{$i[0]}->{e} = $i[10];
        }

        $first{$i[0]} = $i[1] if not defined $first{$i[0]};
        if ($i[1] eq $first{$i[0]} or $i[10] <= $e_thres) {
                push @{$toxin{$i[0]}{al}}, [$i[6], $i[7]];
        }
}
close(FH1);

my %fa;
open FH, $t_file;
my $entry;
my @id_order;
while (<FH>) {
        if (/^>(\S+)/) {
                $entry = $1;
	#debug#	print STDERR "$entry\n";
		push(@id_order, $entry);
        } else {
                chomp;
                $fa{$entry} .= $_;
        }
}
close FH;
my%nPHeaderSeq;
#tie %pepSeq,"Tie::IxHash";
foreach my $tid (@id_order) {
     if(defined $toxin{$tid}){

        $toxin{$tid}->{seq} = $fa{$tid};
        my $length = length $fa{$tid};
        my $peptide = translate_with_frame($toxin{$tid}->{seq}, $toxin{$tid}->{frame});
        $peptide = lc $peptide;
        foreach my $al (@{$toxin{$tid}{al}}) {
                my ($l, $r) = @$al;
                ($l, $r) = ($r,$l) if $l > $r;

                if ($toxin{$tid}{frame} <0) {
                        ($l, $r) = ($length -$r +1, $length - $l + 1);
                }

                my $frame_off = abs($toxin{$tid}->{frame})-1;
                $l = int(($l - $frame_off)/3) +1;
                $r = int(($r - $frame_off)/3 -1)+1;

                $l = 1 if($l < 1);

                my $substr = substr $peptide, $l-1, $r-$l+1;
                $substr = uc $substr;
                substr $peptide, $l-1, $r-$l+1, $substr;
        }

        $/="\n";
	$nPHeaderSeq{$tid}=[$toxin{$tid}{hit},$toxin{$tid}{frame},$toxin{$tid}{seq},$peptide];
     }else{
        $nPHeaderSeq{$tid}=[$fa{$tid}];
    }
}
return \%nPHeaderSeq;
}
#----------------------------------------------------------------------------------------------------------------------------
sub translate_codon {
        my $c = shift;

        $c = uc($c);
        $c =~ s/U/T/g;

        return 'A' if $c eq 'GCT'||$c eq 'GCC'||$c eq 'GCA'||$c eq 'GCG';
        return 'R' if $c eq 'CGT'||$c eq 'CGC'||$c eq 'CGA'||$c eq 'CGG'||$c eq 'AGA'
                ||    $c eq 'AGG';
        return 'N' if $c eq 'AAT'||$c eq 'AAC';
        return 'D' if $c eq 'GAT'||$c eq 'GAC';
        return 'C' if $c eq 'TGT'||$c eq 'TGC';
        return 'Q' if $c eq 'CAA'||$c eq 'CAG';
        return 'E' if $c eq 'GAA'||$c eq 'GAG';
        return 'G' if $c eq 'GGT'||$c eq 'GGC'||$c eq 'GGA'||$c eq 'GGG';
        return 'H' if $c eq 'CAT'||$c eq 'CAC';
        return 'I' if $c eq 'ATT'||$c eq 'ATC'||$c eq 'ATA';
        return 'L' if $c eq 'TTA'||$c eq 'TTG'||$c eq 'CTC'||$c eq 'CTA'||$c eq 'CTG'
                ||    $c eq 'CTT';
        return 'K' if $c eq 'AAA'||$c eq 'AAG';
        return 'M' if $c eq 'ATG';
        return 'F' if $c eq 'TTT'||$c eq 'TTC';
        return 'P' if $c eq 'CCT'||$c eq 'CCC'||$c eq 'CCA'||$c eq 'CCG';
        return 'S' if $c eq 'TCT'||$c eq 'TCC'||$c eq 'TCA'||$c eq 'TCG'||$c eq 'AGT'
                ||    $c eq 'AGC';
        return 'T' if $c eq 'ACT'||$c eq 'ACC'||$c eq 'ACA'||$c eq 'ACG';
        return 'W' if $c eq 'TGG';
        return 'Y' if $c eq 'TAT'||$c eq 'TAC';
        return 'V' if $c eq 'GTT'||$c eq 'GTC'||$c eq 'GTA'||$c eq 'GTG';

        return '*' if $c eq 'TAA'||$c eq 'TGA'||$c eq 'TAG';  # * stands for stop codon.
        return 'O'; # O stands for unknown codon;
}

sub revcom {
        my $seq = shift;

        my @seq = split //, $seq;

        my $revcom_seq = '';
        for (my $i = $#seq; $i >=0; $i--) {
                my $character = $seq[$i];
                $character = 'T' if $character eq 'U';
                $character = 't' if $character eq 'u';

                if ($character eq 'A') {
                        $revcom_seq .= 'T';
                }
                elsif ($character eq 'G') {
                        $revcom_seq .= 'C';
                }
                elsif ($character eq 'C') {
                        $revcom_seq .= 'G';
                }
                elsif ($character eq 'T') {
                        $revcom_seq .= 'A';
                }
                elsif ($character eq 'a') {
                        $revcom_seq .= 't';
                }
                elsif ($character eq 't') {
                        $revcom_seq .= 'a';
                }
                elsif ($character eq 'g') {
                        $revcom_seq .= 'c';
                }
                elsif ($character eq 'c') {
                        $revcom_seq .= 'g';
                }
                elsif ($character eq 'N') {
                        $revcom_seq .= 'N';
                }
                elsif ($character eq 'n') {
                        $revcom_seq .= 'n';
                }
        }
        return $revcom_seq;
}

sub translate_with_frame {
        my $seq = shift;
        my $strand_frame = shift;
	 # Change 9/11/14 to fix strand/frame parsing
	 # my ($strand, $frame) = split //, $strand_frame;
        my ($frame)  = $strand_frame =~ /(\d+)/;
        my ($strand) = $strand_frame =~ /([\+\-])/;
        $strand ||= '+';

        $seq = revcom($seq) if $strand eq '-';
        my @seq = split //, $seq;

        my $peptide = '';
        for (my $i= $frame-1; $i <= $#seq -2; $i+=3) {
                my $codon = $seq[$i].$seq[$i+1].$seq[$i+2];
                my $aa = translate_codon($codon);
                $peptide .= $aa;
        }

        return $peptide;
}

