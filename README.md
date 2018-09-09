# ConusPipe
Here is the pipeline for conotoxin discovery.
The database needed for the pipeline are stored in the archive of server of Yandell lab :/archive02/qli/data/conus_2016/db4ConusPipe/ 
Especially for blastx database, I only put the manually confirmed protein sequence in it. As the manual confirmation follows the new doscovery,
you can espect the blastx datase to grow. 
 to run the main pipeline:
nohup pathToConusPipe/bin/discoveryPipeline.v8.pl pathToConusPipe/doc/sample.config pathToConusPipe/doc/codeCharMwPi
