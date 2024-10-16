files=("g2315_OG0001120" "g2216_OG0001984" "g2173_OG0000896" "g7072_OG0005649")
for file in ${files[@]};do
awk '{if(!/>/) sub(/(TA[GA]|TGA)-*$/," ici &",$0); print}' ./fasta_ortho/${file}.fa | sed 's/ ici [ATGC]\{3\}/nnn/1' > ./fasta_ortho/without_stopcodon/${file}.fa_edited.fa
done