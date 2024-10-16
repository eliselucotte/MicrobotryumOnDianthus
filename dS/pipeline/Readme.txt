PIPELINE TO MAKE THE DS PLOT
#--------------------------#

1) create a folder with the name of your species, encoded in the following format "Mv-cat-C212"

2) copy in this folder the folder the following:
- 0.pipeline.sh
- ./codes/
- macse_v2.05.jar
- yn00_template.ctl

2) modify 0.pipeline.sh, code that launches all the other codes (located in ./codes/) in the right order, with all the dependencies.
You only have to change what is in the first section.
Define: 
- name_species
- path to the folder

3) modify the paths in ./codes/code_macse.sh

3) The R code dS_plot.R can be used as a template to make your plot


Things to be mindful of:
- The code works for a comparison with Mv Lagerheimii. If you want to change the ancestral state to another species, you have to go through the codes
- The first python code (get_ortologs.py) selects the sex chromosomes from Mv lagerheimii, if you want to run it on all chromosome, remove the if condition.
/!\ in the code make_orthologous_fastas.py, the path to the GTF files is encoded: '/home/elise/Assembly/ortholog_reconstruction/0.data/GTF_final/'
/!\ in the code get_orthologs.py, the path to the orthofinder results is encoded: '/home/elise/Assembly/ortholog_reconstruction/1.results/Results_Apr12_6/Orthologues/Orthologues_Mv-lag-1253-A1_final.cds.prot_longest_transcript/'
/!\ in the code make_recap_table_DS.py, the path to the bedfile for Mv-lag is encoded: '/home/elise/Assembly/ortholog_reconstruction/0.data/GTF_final/Mv-lag-1253-A1_final.bed'
You will be able to read these files if you are on server03.


