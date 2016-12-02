#!/usr/bin/perl -w 
use strict;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\nUSAGE: calSpSn4ctrl.pl 5243conotoxin.prob.txt 5243noTox.prob.txt probabilityCutOff total\n";


die($usage) unless $ARGV[2];

my $cono = $ARGV[0];
my $noCono = $ARGV[1];
my $probCutoff = $ARGV[2];
my $total = $ARGV[3];
open(FH,'<',$cono) || die"can't open $cono\n";
my($sp,$sn);
my$noConoFinal=0;
my$conoFinal=0;
while (defined(my $line=<FH>)){
	chomp ($line);
	if($line>=$probCutoff){
		$conoFinal++;
	}
}
close FH;
open(FH,'<',$noCono) || die"can't open $noCono\n";
while (defined(my $line=<FH>)){
	chomp ($line);
	if($line>=$probCutoff){
		$noConoFinal++;
	}
}
close(FH);
$sp=($total-$noConoFinal)/$total;
$sn=$conoFinal/$total;
my$fpr=1-$sp;
my$j=$sp+$sn-1;
print "$probCutoff\t$fpr\t$sn\t$sp\t$j\n";	

	
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
