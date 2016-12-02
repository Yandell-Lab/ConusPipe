#!/usr/bin/perl -w 
use strict;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "USAGE: extractSignalSeqPerSupfam4mafft.pl <cono.gff> <cono.fa> \n fa is in one line header, one line body format";


die($usage) unless $ARGV[1];

my $gffFile = $ARGV[0];
my $faFile = $ARGV[1];
my%sigSeq;
open(FH,'<',$gffFile) || die"can't open $gffFile\n";
while (defined(my $line=<FH>)){
	next if $line=~ /^##/;
		chomp ($line);
		my@h= split /\t/,$line;
		my($supfam)=$h[0]=~ /conotoxin\.(\S+)/;
		$sigSeq{$supfam}{$h[0]}=[$h[3],$h[4]];
}

close(FH);
my(%fa,$cid);
open(FH1,'<',$faFile) || die"can't open $faFile\n";
while (defined(my $line1=<FH1>)){
	chomp($line1);
	if($line1=~ /^>/){
		($cid)=$line1=~ />(\S+)/;

	}else{
		$fa{$cid}=$line1;
	}

}

foreach my $sFam(keys %sigSeq){
	open OUT,">$sFam.fa";
	foreach my$id(keys %{$sigSeq{$sFam}}){
		my$start=$sigSeq{$sFam}{$id}->[0];
		my$end=$sigSeq{$sFam}{$id}->[1];
		my$seq=$fa{$id};
		my$len=$end-$start+1;
		my$signalSeq=substr $seq,$start-1,$len;
		print OUT ">$id\n$signalSeq\n";
		

	}
	close OUT;
}	

	
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
