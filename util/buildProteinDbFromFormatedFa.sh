#! /bin/bash


set -x
FORMATED_FA=$1
KMER_LENGTH=$2
TAXIS_CUT_OFF=$3
KC_OUT_NAME=$4
STI_OUT_NAME=$5
CLASSIFIER_DB_NAME=$6
#SP_SN_OUT=$6

#-z means zero length
if [[ -z $FORMATED_FA || -z $KMER_LENGTH || -z $TAXIS_CUT_OFF || -z $KC_OUT_NAME || -z $STI_OUT_NAME || -z $CLASSIFIER_DB_NAME ]]; then
    echo "Usage: $0 FORMATED_FA	KMER_LENGTH TAXIS_CUT_OFF KC_OUT_NAME STI_OUT_NAME CLASSIFIER_DB_NAME"
    exit 1
#fi closes the if statement
fi
echo -e  "build classifier protein db with given formated fa, kmer length, taxis cut off "
#substitue id of formated protein fa to numbers ->convert protein seq to nt seq ->kanalyze count -> make sti and tri files from protein converted nt seq-> print config file for building db and run build_classification_db 
generateKeySubstId.pl $FORMATED_FA
echo $FORMATED_FA
regex=([A-Za-z0-9._]+)\.fa
#[["$FORMATED_FA" =~ $regex]] && name="${BASH_REMATCH[1]}"
if [[	$FORMATED_FA =~ $regex ]]; then 
 name="${BASH_REMATCH[1]}"
echo "$name"
~/anaconda/bin/python ~/latest_taxonomer/taxonomer/utils/convert_protein_db.py $name.numId.fa $name.p2n.fa
java -jar ~/applica/kanalyze-0.9.7/kanalyze.jar count -k $KMER_LENGTH -d 30 -l 30 -o $KC_OUT_NAME -m hex -f fasta $name.p2n.fa
~/anaconda/bin/python ~/applica/Taxonomer_py/utilities/create_custom_db_2.0.py -i $name.p2n.fa -o $STI_OUT_NAME 

cat >./$name.DB.config <<EOL
db_prefix:$CLASSIFIER_DB_NAME
sti_file:$STI_OUT_NAME.sti
input_fasta:${STI_OUT_NAME}_taxonomer.fa
kmer_length:$KMER_LENGTH
kmer_cutoff:$TAXIS_CUT_OFF
protein:1
afterburner:0
ncpus:15
EOL
~/latest_taxonomer/taxonomer/build_classification_db $name.DB.config
fi
