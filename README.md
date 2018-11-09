# ConusPipe
## Obtain a version of the ConusPipe package from github
ConusPipe is a free pipeline for conotoxin discovery in academic research use.It ensembles 3 machine learning models and a cross-species blastp method to utilize prominent chemical characters of conotoxins to predict whether a certain transcript in a Conus transcriptome, which has no otherwise detectable homologs in current reference databases, is a putative conotoxin.
To download ConusPipe from github, simply open a linux terminal and type:$ git clone https://github.com/Yandell-Lab/ConusPipe
## Usage
In order to run the pipeline, you need to have ncbi blast and signalP installed and in your path. You also need to have python scikit-learn package installed (recommend install it with Anaconda, python 2.7+ version). 
 The full version of pipeline script is ConusPipe/bin/discoveryPipeline.v8.pl. It has certain requirements on the format of input sequence. You can see the example input files (im.16.tri.no_anot.fa and im.16.tri.no_anot.pep.fa) for the full version in ConusPipe/doc.It requires species be labeled before sequence ID of input sequence and tpm value on the header of input nucleotide sequence. The best way to prepare the input files is to run our conotoxin annotation pipeline - tABRA.v2.pl(currently only support paired-end reads) on your raw sequencing reads.
All the input files, databases and output file name should be given in sample.config file (example sample.config file is in ConusPipe/doc) with path. 
To run it, you can copy the sample.config file to your run directory and modify the path and name of the input, output and database files according to your need and type in linux terminal:
nohup pathToConusPipe/ConusPipe/bin/discoveryPipeline.v8.pl sample.config pathToConusPipe/ConusPipe/doc/codeCharMwPi pathToConusPipe/ & 



If you don't want to run our conotoxin annotation pipeline, which starts from fastq files, you can run the simple version of discovery pipeline script ConusPipe/bin/discoveryPipe.4methods.simple.v2.pl, which does not require tpm value on the header and can take your own peptide or nucleotide sequence as input.However, it still requires species be labeled before sequence ID of input sequence and add ".TR" as suffix of species label (you can see the example input file for simple version  ConusPipe/doc/ay.pep.fa).  
To run it, you can type in linux terminal:

nohup pathToConusPipe/ConusPipe/bin/discoveryPipe.4methods.simple.v2.pl [-options] conotoxinseqTraining noConotoxinseqTraining lengthLowCutoff lengthCutoff evalCutoff DvalCutoff  pathToConusPipe/ConusPipe/doc/codeCharMwPi pathToConusPipe/ outName

Options:

--nucleotide <nt.fa> 

--peptide <pep.fa>



The basic trainig database for conotoxinseq and noConotoxinseq are provided in ConusPipe/doc as all.cono.training.pep.fa and conus.notTox.training.50to200Len.pep.fa. As you discover/confirm more new conotoxins and non-conotoxins, you can add more sequence to the training database.  
### For example to run test data with simple version:
nohup pathToConusPipe/ConusPipe/bin/discoveryPipe.4methods.simple.v2.pl --peptide pathToConusPipe/ConusPipe/doc/all.cono.training.pep.fa pathToConusPipe/ConusPipe/doc/conus.notTox.training.50to200Len.pep.fa 6 1000 1e-10 0.45 pathToConusPipe/ConusPipe/doc/codeCharMwPi pathToConusPipe/ testSimple & 
## Expected output:
For each method, it outputs the overlapped sequence in xx.(4,3,2).overlap.pep.fa.txt, print out non-overlapped sequence in xx.bp.pep.fa.txt, xx.logit.pep.fa.txt, xx.labelSpread.pep.fa.txt, xx.pep.cono.fa.txt. Finally all the predicted putative toxins are put together in the file xx.total.ml.pep.fa. 
For the full version pipeline, there's also a tblastn run against ncbi-nr molluscan database (2018/03 version), and if there's a hit, the result is in xx.mtbn.nt.fa and xx.mtbn.pep.fa. 





