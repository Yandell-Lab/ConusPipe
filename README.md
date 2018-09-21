# ConusPipe
Here is the pipeline for conotoxin discovery.
The database and training data needed for the pipeline are stored in the archive of server of Yandell lab :/archive02/qli/data/conus_2016/db4ConusPipe/ 
Especially for blastx database, I only put the manually confirmed protein sequence in it. As the manual confirmation follows the new doscovery,
you can espect the blastx datase to grow. 
 to run the main pipeline (the input peptide and nucleotide sequence should be output from our conotoxin annotation pipeline - tABRA.v2.pl) to fulfill certain format requirement:
nohup pathToConusPipe/bin/discoveryPipeline.v8.pl pathToConusPipe/doc/sample.config pathToConusPipe/doc/codeCharMwPi

If you don't want to run our conotoxin annotation pipeline, which starts from fastq files, you can run the simple version of discovery pipeline, which does not have certain format requirement and can take your own peptide or nucleotide sequence as input:
nohup pathToConusPipe/bin/discoveryPipe.4methods.simple.pl [-options] conotoxinseqTraining noConotoxinseqTraining fpDb  lengthLowCutoff lengthCutoff evalCutoff DvalCutoff  pathToConusPipe/doc/codeCharMwPi outName
Options:
--nucleotide <nt.fa>
--peptide <pep.fa>



