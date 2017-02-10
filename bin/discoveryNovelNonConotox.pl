#!/usr/bin/perl
use strict;
use warnings;
use Tie::IxHash;
use Bio::DB::Fasta;
use Statistics::R;
use Getopt::Long;
# coded by Qing Li on 07/29/16.
#Copyright (c) 2016 yandell lab. All rights reserved.


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "discoveryNovelNonConotox.pl  xx.transcriptome.nt.fa xx.bx.out  tpmCutoff outName\n
input xx.transcriptome.nt.fa xx.bx.out should be the output of previous conotoxin annotation pipeline - xx.transcriptome.nt.fa should have tpm value on the header.
#1.read the xx.transcriptome.nt.fa, skip every entry with tpm < tpmCutoff, store the header and nt seq in hash of array ref;{id}->[Entry,len,tpm,seq]
#2. for transcripts with no hits: 6 frame translation of assembly output,extract ORF >50aa; for transcripts with hits, extract ORF(>50aa) accoridng to frame information in bx.out
#3. for all ORFs, remove the redundant seq(only keep the longest),since signalp can only start to predict from the 1st M, so have to trim the seq according to different Ms. for the trimmed ones,filtering by requiring length > 50aa. Then run signalP( -f short) on pepseq to get output, parse output(D-value > 0.7,9th column==Y) to extract potential new toxins from pepORF,remove redundant ones(only keep the longest),trim the nt.ORF.fa accordingly, put into hash {id}->[ntOrf,pepOrf];
#4. sort ntHash{id}->[Entry,len,tpm,seq] according to tpm value
#5. output i) the extracted complete nt sequences in a fasta file, (ii) the extracted ORF nt sequences in a fasta file and (iii) the extracted ORF translated peptide sequences in a fasta file according to the sorted order of ntHash. 
";
my$ntFa=$ARGV[0];
my$bxOut=$ARGV[1];
my $tpmCutoff=$ARGV[2];
my $sample=$ARGV[3];
die $usage unless $sample;
my%frame;
open FH1,$bxOut;
while(defined (my$line=<FH1>)){
        chomp ($line);
	my@i=split /\t/,$line;
	$frame{$i[0]}=$i[12];
}
close FH1;
my(%ntHash,$id,$entry,$len,$tpm,$seq);#{id}->[Entry,len,tpm,seq] 
my%pepOrfHash;#{pepOrf}->$id
open OUT,">$sample.complete.pep.fa"; 
open FH,$ntFa; 
while(defined (my$line=<FH>)){
	chomp ($line);
	if($line=~ />/){
        	my@t=split /\t/,$line;
		$t[0]=~ s/>//;
		$t[1]=~ s/Entry://;
		$t[2]=~ s/len://;
		$t[3]=~ s/tpm://;
		$id=$t[0];
		$entry=$t[1];
		$len=$t[2];
		$tpm=$t[3];
	}else{
		next if $tpm < $tpmCutoff;
		$seq=$line;	
		$ntHash{$id}=[$entry,$len,$tpm,$seq];
		if($entry=~ /noBxHit/){
			my@frames=(-1,-2,-3,1,2,3);
			 for my $j (@frames) {
	                        translateSeq($id,$seq, $j,\%pepOrfHash,\*OUT);
			}
                }else{
			my$fr=$frame{$id};
			translateSeq($id,$seq, $fr,\%pepOrfHash,\*OUT);
		}
		
	}
}
close FH;
close OUT;
## for all ORFs, remove the redundant seq(only keep the longest).
my@pepOrf= keys %pepOrfHash;
my$uniqOrf=removeRedundant(\@pepOrf);
#since signalp can only start to predict from the 1st M, so have to trim the seq according to different Ms. for the trimmed ones,filtering by requiring length > 50aa.
my %sigPNumId; #{numId}->seqId with trimming info
my$idCount=1;
open OUT8,">$sample.beforeSigP.pep.fa";
foreach my$uSeq(@$uniqOrf){
	my$offset=0;
	my$mIndex=index($uSeq,'M',$offset);	
	my$count=0;
	while($mIndex !=-1){
		my$uSeqM=substr $uSeq,$mIndex;	
		my$checkLen=length $uSeqM;
		last unless $checkLen > 50 ;
		my$idLen=length("$pepOrfHash{$uSeq}.trim$count");
		if($idLen<55){
			print OUT8 ">$pepOrfHash{$uSeq}.trim$count\n$uSeqM\n" ;
		}else{
			print OUT8 ">$idCount\n$uSeqM\n";
                        $sigPNumId{$idCount}="$pepOrfHash{$uSeq}.trim$count";
		}
		$offset=$mIndex+1;
		$count++;
		$mIndex=index($uSeq,'M',$offset);
		$idCount++;
	}
}
close OUT8;
#run signalp
system("signalp -t euk -u 0.7 -U 0.7 -f short  $sample.beforeSigP.pep.fa > $sample.sigP.short_out");
#parse signalp output 
my%sigPSeq;
my($beforeSigPFa)=<$sample.beforeSigP.pep.fa>;
my$SigPDb=Bio::DB::Fasta->new($beforeSigPFa);
open(FH2, "$sample.sigP.short_out")or die "Can't open $sample.sigP.short_out for reading: $!\n";
while(defined(my$line=<FH2>)){
	next if $line=~ /^#/;
	my@p=split /\s+/,$line;
	next if $p[9] eq "N";
	my$seq=$SigPDb-> get_Seq_by_id($p[0]);
        my $seqstr  = $seq-> seq;
	if(!defined $sigPNumId{$p[0]}){
                $sigPSeq{$seqstr}=$p[0];
        }else{
                $sigPSeq{$seqstr}=$sigPNumId{$p[0]};
        }

}
close FH2;
#remove redundant ones(only keep the longest),trim the nt.ORF.fa accordingly, put into hash {id}->[ntOrf,pepOrf];
my%sigPIdSeq;#{id}->[tpm,ntOrf,pepOrf]
my@sigPSeq=keys %sigPSeq;
my$uniqSigPSeq=removeRedundant(\@sigPSeq);
my($compPepFa)=<$sample.complete.pep.fa>;
my$nAPDb=Bio::DB::Fasta->new($compPepFa);
foreach my$uSigPSeq(@$uniqSigPSeq){
		my$len=length($uSigPSeq);
		my($dbNId)=$sigPSeq{$uSigPSeq}=~ /(\S+)_frame\S+/;
		my($dbPId)=$sigPSeq{$uSigPSeq}=~ /(\S+)\.trim\S+/;
		my$ntSeq=$ntHash{$dbNId}->[3];
		#print "$ntSeq\n";		
		my$pepSeq=$nAPDb-> get_Seq_by_id($dbPId);
                my $pepSeqstr  = $pepSeq-> seq;
		my ($frame_strand) = $sigPSeq{$uSigPSeq}=~ /frame_(\S+)_/;
                my($frame)  = $frame_strand =~ /(\d+)/;
                my($strand) = $frame_strand =~ /([\+\-])/;
                $strand ||= '+';
		my $mIndex=index($pepSeqstr,$uSigPSeq);
                my$pepLen=length($uSigPSeq);
                my $ntStart=$mIndex*3+$frame-1;
                my $ntLen=$pepLen*3;
		$ntSeq= revcom($ntSeq) if $strand eq '-';  
		my $ntOrf=substr $ntSeq, $ntStart, $ntLen;
		my$tpm=$ntHash{$dbNId}->[2];
		$sigPIdSeq{$sigPSeq{$uSigPSeq}}=[$tpm,$ntOrf,$uSigPSeq];
}
#sort sigPIdSeq{id}->[$tpm,$ntOrf,$uSigPSeq] according to tpm value, output i) the extracted complete nt sequences in a fasta file, (ii) the extracted ORF nt sequences in a fasta file and (iii) the extracted ORF translated peptide sequences in a fasta file according to the sorted order, with the header from ntHash.
my%nidFlag;
my@ids=sort{$sigPIdSeq{$b}->[0] <=> $sigPIdSeq{$a}->[0]} keys(%sigPIdSeq);
open OUT2,">$sample.extract.complete.nt.fa";
open OUT3,">$sample.extract.orf.nt.fa";
open OUT4,">$sample.extract.orf.pep.fa";
foreach my $id (@ids) {
	my($nId)=$id=~ /(\S+)_frame\S+/;	
	my$oriHeader=">$nId\tEntry:$ntHash{$nId}->[0]\tlen:$ntHash{$nId}->[1]\ttpm:$ntHash{$nId}->[2]";
	my$header=">$id\tEntry:$ntHash{$nId}->[0]\tlen:$ntHash{$nId}->[1]\ttpm:$ntHash{$nId}->[2]";
	my$ntSeq=$ntHash{$nId}->[3];
	my$ntOrf=$sigPIdSeq{$id}->[1];
	my$pepOrf=$sigPIdSeq{$id}->[2];
	print OUT3 "$header\n$ntOrf\n";
	print OUT4 "$header\n$pepOrf\n";
	next if defined $nidFlag{$nId}; 
	print OUT2 "$oriHeader\n$ntSeq\n";
	$nidFlag{$nId}=1;
}
close OUT2;
close OUT3;
close OUT4;



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
# translateSeq($id,$seq, $fr,$pepOrfHash,\*OUT)
sub translateSeq{
	my$id=shift;
	my$seq=shift;
	my$frame=shift;
	my$pepOrfHash=shift;
	 my $fh=shift;
	my $pep = translate_with_frame($seq, $frame);
	my @peps=$pep=~ /([A-Z]{50,})\*/g;
                for (my$i=0;$i<@peps;$i++){
			$pepOrfHash->{$peps[$i]}=$id."_frame_$frame"."_$i";
			print $fh ">$id"."_frame_$frame"."_$i\n$peps[$i]\n";
		}
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
#	my ($anot_nt_fa) = <$anot_fa_file.nt>;
#	return $anot_nt_fa;
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

