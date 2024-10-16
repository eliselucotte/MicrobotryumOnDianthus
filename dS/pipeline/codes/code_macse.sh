#!/bin/bash
ulimit -s 1000000
filename=$1
#echo $filename
path='/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/Mv-cat-C212/fasta_ortho/without_stopcodon/'
path_res='/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/Mv-cat-C212/alignement/'

if test -f "${path_res}${filename}_NT.fa"; then #if the file exists, do nothing
echo "${filename}_NT.fa exists"
else
java -jar ./macse_v2.05.jar -prog alignSequences -seq ${path}${filename} -out_NT ${path_res}${filename}_NT.fa -out_AA ${path_res}${filename}_AA.fa 2>&1 | tee ${path_res}/std/${filename}.log
fi