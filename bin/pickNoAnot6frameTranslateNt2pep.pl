#!/usr/bin/perl
use strict;
use warnings;
use Tie::IxHash;
use Bio::DB::Fasta;
use Array::Utils qw(:all);
use Getopt::Long;
use File::Basename;
# coded by Qing Li on 07/29/16.
#Copyright (c) 2016 yandell lab. All rights reserved.


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "pickNoAnot6frameTranslate.pl xx.transcriptome.nt.fa lengthLowCutoff outName\n";
my($triFa,$lengthLowCutoff,$sample)=@ARGV;
die $usage unless $sample;
##6 frame translation
	my$db=Bio::DB::Fasta->new($triFa);
	open FH,$triFa; 
	open OUT6,">$sample.tri.no_anot.fa";
	open OUT7,">$sample.tri.no_anot.pep.fa";
	while(defined (my$line=<FH>)){
		if($line=~ />/){
			my($id)=$line=~ />(\S+)/;
			my$desc=$db->header($id);
			next unless $desc=~ /Entry:noBxHit/;
			my$seqstr=$db->seq($id);
			print OUT6 ">$desc\n$seqstr\n";
			for my $j (1..3){
                        	translate_print_seq($id,$seqstr, $j,\*OUT7,$lengthLowCutoff);
                	}

		}
	}
	close FH;
	close OUT6;
	close OUT7;
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
sub numId{
        my$input=shift;
        open FH1,$input;
	my ($fileName)=basename($input);
        open OUT1,">numId.$fileName";
        my$idCount=1;
        while (defined (my $line=<FH1>)){
                if($line=~ />/){
                        my($id)=$line=~ />(\S+)/;
                        my$idLen=length($id);
                        if($idLen<55){
                                print OUT1 $line;
                        }else{
                                print OUT1 ">$idCount\n";
                                $idCount++;
                        }
                }else{
                        print OUT1 $line;
                }
        }
close OUT1;
}


sub getFeatures{
        my$seq=shift;
        my$sigpOut=shift;
        my$aaInfo=shift;
        my$sigDVal=getDVal($sigpOut);
        open FH,$seq;
        my($id,%features);#{id}->[cys%(ss),cys%(pp),cys%(mtox),molecular weight(ss),molecular weight(pp),molecular weight(mtox),+%(ss),+%(pp),+%(mtox),
# -%(ss),-%(pp),-%(mtox),pi(ss),pi(pp),pi(mtox)]
	while(defined(my$line=<FH>)){
                chomp $line;
                if($line=~ />/){
                        ($id)=$line=~ />(\S+)/;

                }else{
			my($cysSs,$cysPp,$cysMTox);
                        my$len=length($line);
                        my@aa=split //,$line;
                        my($countSs,$countPp,$countMTox,$mwSs,$mwPp,$mwMTox,$pcSs,$pcPp,$pcMTox,$ncSs,$ncPp,$ncMTox,$piSs,$piPp,$piMTox)=0;
                        for(my$i=0;$i<$len-1;$i++){
                                if($i>=0 && $i<$len*0.31-1){
                                        $mwSs+=$aaInfo->{$aa[$i]}->[1];
                                        $piSs+=$aaInfo->{$aa[$i]}->[2];
					$cysSs++ if $aa[$i] eq "C";
                                        $countSs++;
                                        if($aaInfo->{$aa[$i]}->[0] eq "positive"){
                                                $pcSs++;
                                        }elsif($aaInfo->{$aa[$i]}->[0] eq "negative"){
                                                $ncSs++;
                                        }
                                }elsif($i>=$len*0.31-1 && $i<$len*0.58-1){
                                        $mwPp+=$aaInfo->{$aa[$i]}->[1];
                                        $piPp+=$aaInfo->{$aa[$i]}->[2];
					$cysPp++ if $aa[$i] eq "C";
                                        $countPp++;
                                        if($aaInfo->{$aa[$i]}->[0] eq "positive"){
                                                $pcPp++;
                                        }elsif($aaInfo->{$aa[$i]}->[0] eq "negative"){
                                                $ncPp++;
                                        }

                                }else{
                                        $mwMTox+=$aaInfo->{$aa[$i]}->[1];
                                        $piMTox+=$aaInfo->{$aa[$i]}->[2];
					$cysMTox++ if $aa[$i] eq "C";
                                        $countMTox++;
                                        if($aaInfo->{$aa[$i]}->[0] eq "positive"){
                                                $pcMTox++;
                                        }elsif($aaInfo->{$aa[$i]}->[0] eq "negative"){
                                                $ncMTox++;
                                        }
                                }

                        }
			my$cysSsPct=$cysSs/$countSs;
			my$mwSsAve=$mwSs/($countSs*205);
                        my$piSsAve=$piSs/($countSs*11);
                        my$pcSsPct=$pcSs/$countSs;
                        my$ncSsPct=$ncSs/$countSs;
			my$cysPpPct=$cysPp/$countPp;
                        my$mwPpAve=$mwPp/($countPp*205);
                        my$piPpAve=$piPp/($countPp*11);
                        my$pcPpPct=$pcPp/$countPp;
                        my$ncPpPct=$ncPp/$countPp;
			my$cysMToxPct=$cysMTox/$countMTox;
                        my$mwMToxAve=$mwMTox/($countMTox*205);
                        my$piMToxAve=$piMTox/($countMTox*11);
                        my$pcMToxPct=$pcMTox/$countMTox;
                        my$ncMToxPct=$ncMTox/$countMTox;
                        if($sigDVal->{$id}){
                                $features{$id}=[$cysSsPct,$cysPpPct,$cysMToxPct,$sigDVal->{$id},$mwSsAve,$mwPpAve,$mwMToxAve,$pcSsPct,$pcPpPct,$pcMToxPct,$ncSsPct,$ncPpPct,$ncMToxPct,$piSsAve,$piPpAve,$piMToxAve];
          		}else{
                                $features{$id}=[$cysSsPct,$cysPpPct,$cysMToxPct,0,$mwSsAve,$mwPpAve,$mwMToxAve,$pcSsPct,$pcPpPct,$pcMToxPct,$ncSsPct,$ncPpPct,$ncMToxPct,$piSsAve,$piPpAve,$piMToxAve];
			}
		}

	}
	return \%features;
        close FH;
}






sub getXFeatures{
	my$sigPIdSeq=shift;
        my$aaInfo=shift;
        my$sigDVal=shift;
	my%features;
	foreach my $sid (keys %{$sigPIdSeq}){
        	if(1>0){#change the structure of bp to parallel
			#my@cys=$sigPIdSeq->{$sid}->[3]=~ /(C)/g;
			my($cysSs,$cysPp,$cysMTox);
			my$len=length($sigPIdSeq->{$sid});
			my@aa=split //, $sigPIdSeq->{$sid};
			my($countSs,$countPp,$countMTox,$mwSs,$mwPp,$mwMTox,$pcSs,$pcPp,$pcMTox,$ncSs,$ncPp,$ncMTox,$piSs,$piPp,$piMTox);
                        for(my$i=0;$i<$len-1;$i++){
                                if($i>=0 && $i<$len*0.31-1){
                                        $mwSs+=$aaInfo->{$aa[$i]}->[1];
                                        $piSs+=$aaInfo->{$aa[$i]}->[2];
					$cysSs++ if $aa[$i] eq "C";
                                        $countSs++;
                                        if($aaInfo->{$aa[$i]}->[0] eq "positive"){
                                                $pcSs++;
                                        }elsif($aaInfo->{$aa[$i]}->[0] eq "negative"){
                                                $ncSs++;
                                        }
                                }elsif($i>=$len*0.31-1 && $i<$len*0.58-1){
                                        $mwPp+=$aaInfo->{$aa[$i]}->[1];
                                        $piPp+=$aaInfo->{$aa[$i]}->[2];
					$cysPp++ if $aa[$i] eq "C";
                                        $countPp++;
                                        if($aaInfo->{$aa[$i]}->[0] eq "positive"){
                                                $pcPp++;
                                        }elsif($aaInfo->{$aa[$i]}->[0] eq "negative"){
                                                $ncPp++;
                                        }

                                }else{
                                        $mwMTox+=$aaInfo->{$aa[$i]}->[1];
                                        $piMTox+=$aaInfo->{$aa[$i]}->[2];
					$cysMTox++ if $aa[$i] eq "C";
                                        $countMTox++;
                                        if($aaInfo->{$aa[$i]}->[0] eq "positive"){
                                                $pcMTox++;
                                        }elsif($aaInfo->{$aa[$i]}->[0] eq "negative"){
                                                $ncMTox++;
                                        }
                                }

                        }
                        #print "$mwSs\t$piSs\t$pcSs\t$ncSs\t$countSs\t$mwPp\t$piPp\t$pcPp\t$ncPp\t$countPp\t$mwMTox\t$piMTox\t$pcMTox\t$ncMTox\t$countMTox\n";#debugger
                        my$cysSsPct=$cysSs/$countSs;
                        my$mwSsAve=$mwSs/($countSs*205);
                        my$piSsAve=$piSs/($countSs*11);
                        my$pcSsPct=$pcSs/$countSs;
                        my$ncSsPct=$ncSs/$countSs;
			my$cysPpPct=$cysPp/$countPp;
                        my$mwPpAve=$mwPp/($countPp*205);
                        my$piPpAve=$piPp/($countPp*11);
                        my$pcPpPct=$pcPp/$countPp;
                        my$ncPpPct=$ncPp/$countPp;
			my$cysMToxPct=$cysMTox/$countMTox;
                        my$mwMToxAve=$mwMTox/($countMTox*205);
                        my$piMToxAve=$piMTox/($countMTox*11);
                        my$pcMToxPct=$pcMTox/$countMTox;
                        my$ncMToxPct=$ncMTox/$countMTox;
			if($sigDVal->{$sid}){
                                $features{$sid}=[$cysSsPct,$cysPpPct,$cysMToxPct,$sigDVal->{$sid},$mwSsAve,$mwPpAve,$mwMToxAve,$pcSsPct,$pcPpPct,$pcMToxPct,$ncSsPct,$ncPpPct,$ncMToxPct,$piSsAve,$piPpAve,$piMToxAve];
			}else{
                                $features{$sid}=[$cysSsPct,$cysPpPct,$cysMToxPct,0,$mwSsAve,$mwPpAve,$mwMToxAve,$pcSsPct,$pcPpPct,$pcMToxPct,$ncSsPct,$ncPpPct,$ncMToxPct,$piSsAve,$piPpAve,$piMToxAve];
			}
		}
	}
	return \%features;

}
sub getHit{
	my $blast_file = shift;
	open FH1,$blast_file;
	my %toxin;
	while (<FH1>) {
		my @i = split /\t/, $_;
		if (not defined $toxin{$i[0]}) {
			$toxin{$i[0]} = $i[1];
		}
	}
	return (\%toxin);	

}
sub cleanSpName{
        my$sp=shift;
        $sp=~ s/SG//;
        $sp=~ s/VD//;
        $sp=~ s/NR//;
        $sp=~ s/_MP//;
        $sp=~ s/4section//;
        return $sp;
}


sub getCysDL{
	my$seq=shift;
	my$sigpOut=shift;
	my$sigDVal=getDVal($sigpOut);
	open FH,$seq;
	my($id,%cysD);#%cysD:{id}->[cys%, D-val,len]
	while(defined(my$line=<FH>)){
		chomp $line;
		if($line=~ />/){
			($id)=$line=~ />(\S+)/;

		}else{
			my@cys=$line=~ /(C)/g;
                	my$len=length($line);
			my$lenNorm=$len/784;
                	my$cysNum=scalar @cys;
                	my$cysPct=$cysNum/$len;	
			if($sigDVal->{$id}){
				$cysD{$id}=[$cysPct,$sigDVal->{$id},$lenNorm];
			}else{
				$cysD{$id}=[$cysPct,0,$lenNorm];
			}
		}	


	}
	
	return \%cysD;
	close FH;
}

sub getDVal{
	my $sigpOut=shift;
	my%sigDVal;
	open(FH2, $sigpOut)or die "Can't open $sigpOut for reading: $!\n";
	while(defined(my$line=<FH2>)){
		next if $line=~ /^#/;
		my@p=split /\s+/,$line;
		$sigDVal{$p[0]}=$p[8]; 
	}
	return \%sigDVal;
	close FH2;
}
sub calculateN50{
	my$sample=shift;
	my($filename)=<$sample.trinity.Trinity.fasta>;
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






sub translate_print_seq {
        my $entry=shift;
        my $seq=shift;
        my $frame=shift;
	my $fh=shift;
	my $pep_RL=shift;
        my $pep = translate_with_frame($seq, $frame);
        my $rpep = translate_with_frame($seq, -$frame);
	#print $fh $entry."_frame_$frame\n$pep\n";
	#print $fh $entry."_frame_-$frame\n$rpep\n";

        my @peps=$pep=~ /([A-Z]{$pep_RL,})\*/g;
                for (my$i=0;$i<@peps;$i++){
                        print $fh ">$entry"."_frame_$frame"."_$i\n$peps[$i]\n";
                }

        my @rpeps=$rpep=~ /([A-Z]{$pep_RL,})\*/g;
                for (my$k=0;$k<@rpeps;$k++){
                        print $fh ">$entry"."_frame_-$frame"."_$k\n$rpeps[$k]\n";
                }

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
	close OUT1;
	close FH1;

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
my $order=1;
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

