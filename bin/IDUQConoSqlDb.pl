#!/usr/bin/perl
use DBI;
use Getopt::Long;
use strict;
my$usage="IDUQConoSqlDb.pl [-options] 
Options:
--insert <newCono.fa>
--delete <seqId.txt>
--update <seqId.txt>
--querySF <supFam.txt>
--querySP <species.txt>
--db	<database>
note: fasta file should be one line header, one line body format
       the *.txt file should be single string/line
	database file is required!\n";
my($newConoFa,$seqIdFile4d,$seqIdFile4u,$supFamFile,$spFile,$conoDb);
GetOptions('insert=s'	=>\$newConoFa,
	   'delete=s'   =>\$seqIdFile4d,
	   'update=s'	=>\$seqIdFile4u,	
	   'querySF=s'    =>\$supFamFile,
           'querySP=s'    =>\$spFile, 
	   'db=s'	  =>\$conoDb
);
die $usage unless $conoDb;
my $driver   = "SQLite"; 
my $dsn = "DBI:$driver:dbname=$conoDb";
my $userid = "";
my $password = "";
my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 }) 
                      or die $DBI::errstr;

print "Opened database successfully\n";
insert($newConoFa,$dbh) if $newConoFa;
deletes($seqIdFile4d,$dbh) if $seqIdFile4d;
update($seqIdFile4u,$dbh) if $seqIdFile4u;
querySF($supFamFile,$dbh) if $supFamFile;
querySP($spFile,$dbh) if $spFile;
#-------------------------------------------------SUBS------------------------------------------------------------------------------------------------------------------------------------------

sub insert {
	my$newConoFa=shift;
	my$dbh=shift;
	my($stmt,%idHash);
	my($id,$spInfo,$supfam);
	open FH1,$newConoFa;
	while(defined (my$line=<FH1>)){	
		chomp $line;
		my$seq;
		if($line=~ />/){
			($id)=$line=~ />(\S+)/;
			if($id=~ /TRINITY/){
				($spInfo)=$id=~ /(\S+)\.TRINITY/;
				my@header=split /\t/,$line;
				$header[1]=~ s/supFam://;
				$header[1]=~ s/putative//;
				$supfam=$header[1];
			}else{
				die "fasta file is not in the right format!!\n"
			}
		}else{
			$seq=$line;
		}
		if ($seq=~ /\w+/){
			next if defined $idHash{$id};
			$idHash{$id}=1;
			$stmt = $dbh->prepare('INSERT INTO conoSupfamSp VALUES (?,?,?,?)');
			$stmt->execute($id,$supfam,$seq,$spInfo) or die $DBI::errstr;
			print "Records created successfully\n";
		}
	}
	$dbh->disconnect();
}
	
sub deletes{
	my $seqIdFile=shift;	
	my$dbh=shift;	
	open FH1,$seqIdFile;
	while(defined (my$id=<FH1>)){
		chomp $id;
		my$stmt=$dbh->prepare('DELETE from conoSupfamSp where ID = ?');
		$stmt->execute($id) or die $DBI::errstr;
		print "Records deleted successfully\n"
	}
	$dbh->disconnect();	
}

sub update{
	my $seqIdFile=shift;	
	my$dbh=shift;	
	open FH1,$seqIdFile;
	while(defined (my$id=<FH1>)){
		chomp $id;
		my$stmt=$dbh->prepare('UPDATE conoSupfamSp set SPINFOR= new where ID = ?');
		$stmt->execute($id) or die $DBI::errstr;
		print "Records updated successfully\n"
	}
	$dbh->disconnect();	
}
sub querySF{
	my $supFam=shift;
        my$dbh=shift;
        open FH1,$supFam;
        print "id\tsupfam\tseq\tspInfo\n";
        while(defined (my$sf=<FH1>)){
                chomp $sf;
                my$stmt=$dbh->prepare('SELECT * from conoSupfamSp where SUPFAM = ? COLLATE NOCASE;');
                $stmt->execute($sf) or die $DBI::errstr;
		while(my @row = $stmt->fetchrow_array()) {
			print"$row[0]\t$row[1]\t$row[2]\t$row[3]\n";
		}
        }
        $dbh->disconnect();
}

sub querySP{
	my $spFile=shift;
        my$dbh=shift;
        open FH1,$spFile;
        print "id\tsupfam\tseq\tspInfo\n";
        while(defined (my$sp=<FH1>)){
                chomp $sp;
                my$stmt=$dbh->prepare('SELECT * from conoSupfamSp where SPINFOR like ? COLLATE NOCASE;');
                $stmt->execute($sp.'%') or die $DBI::errstr;
		while(my @row = $stmt->fetchrow_array()) {
			print"$row[0]\t$row[1]\t$row[2]\t$row[3]\n";
		}
        }
        $dbh->disconnect();
}
