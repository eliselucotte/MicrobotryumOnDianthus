###############################
# PIPELINE TO CALCULATE THE DS#
###############################
## Author : Elise Lucotte
## date of last modification : June 2023
## description : launches the different step to calculate the DS between 
## Mv Lagerheimii and another species
###############################

## Launch the pipeline and keep a log:
## ./0.pipeline.sh 2>&1 |tee log.pipeline
name_species='Mv-sup-1065'
##Path to the working directory
path=/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/${name_species}/
##Path to the orthofinder results
path_ortho=/home/elise/Assembly/ortholog_reconstruction/1.results/Results_Jun12_2/
##0) Creates the folder structure
##--------------------------------
mkdir ${path}/fasta_ortho/
mkdir ${path}/fasta_ortho/without_stopcodon/
mkdir ${path}/alignement/
mkdir ${path}/alignement/std/
mkdir ${path}/alignement/correspondance_id/
mkdir ${path}/paml/
mkdir ${path}/paml/results/
mv yn00_template.ctl ${path}/paml/

##1) Prepare the orthologs fastas
##-------------------------------
echo 'prepare the ortologs fasta'
##Creates a files with the info of the ortologs of Mvlag on the sex chromosomes, for the species of interest
python3 ${path}/codes/get_orthologs.py $name_species $path_ortho

##Creates the fasta for each ortologous pair A1 A2 for the species of interest
python3 ${path}/codes/make_orthologous_fastas.py $name_species

## Removes the stop codon from the fasta files
cd ${path}/fasta_ortho/
for file in $(ls *.fa); do
	#if [ -s $file ]; then echo $file; fi #echo the not empty files
	cat ${path}/fasta_ortho/${file} | awk '{if($1 !~ /^>/) { {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)}} } else{print $1}}' > ${path}/fasta_ortho/without_stopcodon/${file}_edited.fa
done

##2) Alignement#
##-------------#
echo 'Alignments'
##make a list of files sorted by size
ls -Sr ${path}/fasta_ortho/without_stopcodon/ > ${path}list_files_ordered.txt
cd ${path}
##launches macse 20 jobs in parallel
##has a function inside that looks at whether the output file already exists
parallel --jobs 20 -a ${path}list_files_ordered.txt ${path}/codes/code_macse.sh $1

## Reformating the alignments for PAML
cd ${path}/alignement/
for file in $(ls *NT.fa); do
	ortho=${file%%.*}
	echo $ortho
	##check if the reformating was already done
	final_file=${ortho}_NT_phylip.fa
	if test -f "$final_file"; then #if the file exists, do nothing
    echo "$final_file exists"
	else
	##remove the ! from the alignments
	echo 'removing the !'
	sed 's/!/-/g' ${ortho}.fa_edited.fa_NT.fa > ${ortho}_NT_noexcl.fa

	##change the ids into number and create a table with the correspondance into the folder correspondance_id/
	echo 'changing the ID'
	python3 ${path}/codes/change_ids.py $path $ortho 

	##convert the alignment into a phylip format
	echo 'convert into phylip format'
	python3 ${path}/codes/convert_macse_to_phylip.py $path $ortho

	##interleave the phylip file
	echo 'interleaving the phylip file'
	sed -i '1 s/$/ I/' ${path}/alignement/${ortho}_NT_phylip.fa
	fi
done

##3- PAML #
##--------#
echo 'PAML'
cd ${path}/alignement/
path_paml=${path}/paml/
path_alignment=${path}/alignement/
##Making the control files#

for file in $(ls *NT_phylip.fa); do
	##extract the name of the ortholog from the file name
	ortho=${file%%_NT*}
	##check if the final reformated file already exists
	ctl_file=${path_paml}/results/${ortho}/yn00_${ortho}.ctl
	if test -f "$ctl_file"; then
    echo "$ctl_file exists"
	
	else
	echo $ortho
	##else, create the folder and run the python script to create the ctl file
	mkdir ${path_paml}/results/${ortho}/
	python3 ${path}/codes/make_control_files.py ${path_alignment} ${path_paml} ${ortho}
	fi
done

##launch yn00
##-----------
cd ${path}/alignement/
for file in $(ls *NT_phylip.fa); do
	##extract the name of the ortholog from the file name
	ortho=${file%%_NT*}
	##check if the final file already exists
	outfile=../paml/results/${ortho}/yn00_orthogp.${ortho}.out
	if test -f "$outfile"; then
    echo "$outfile exists"
	
	else
	echo $ortho
	##else, run the ctl file 
	yn00 ../paml/results/${ortho}/yn00_${ortho}.ctl
	mv ${path}/paml/yn00_orthogp.${ortho}.out ${path}/paml/results/${ortho}/
	fi
done


##Parse the results
##-----------------
python3 ${path}/codes/parse_results_paml_DS.py ${path}
python3 ${path}/codes/make_recap_table_DS.py ${path} ${name_species}



