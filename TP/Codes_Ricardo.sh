for i in $(ls -1 /home/elise/Assembly/transpecific_polymorphism/trees/*/*iqtree \
| awk -F "/" '{print $NF}'); do awk '/;$/{split(FILENAME,f,"/"); gsub(/_NT.+$/,"",f[8]); print f[8],f[7],$1}' \
/home/elise/Assembly/transpecific_polymorphism/trees/*/$i; done \
| sed -f <(ls /home/elise/Assembly/transpecific_polymorphism/alignement/correspondance_id/*txt \
| awk -F "/" '{split($NF,f,"_"); print "awk xx{print \"/"f[1]"_"f[2]"/ s/\"$NF\":/\"$2\":/1\"}xx "$0}' \
| sed "s/xx/'/g" | bash) - > g2ogf2stratum2tree

/home/ricardo/Tools/newick_utils/src/nw_order <(awk '{print $NF}' g2ogf2stratum2tree) \
| awk '{printf ">%03i %s\n",NR,$0}' \
| sort -k1,1 \
| join -j1 <(awk '{printf ">%03i %s\n",NR,$0}' g2ogf2stratum2tree \
	| awk '{print $1,$2,$3}' | sort -k1,1) - \
| sed 's/_h2tg000005l_g[0-9]\+.t[0-9]//g;s/_MC03_g[0-9]\+.t[0-9]//g;s/_tig00000080_g[0-9]\+.t[0-9]//g;s/_h1tg000006l_g[0-9]\+.t[0-9]//g;s/_tig00000068_g[0-9]\+.t[0-9]//g;s/_tig00000084_g[0-9]\+.t[0-9]//g' \
> g2ogf2stratum2tree_reordered

awk '{if(/Mv-cat-C212-A1:[0-9]+\.[0-9]+,Mv-sup-1065-A1/ && /Mv-cat-C212-A2:[0-9]+\.[0-9]+,Mv-sup-1065-A2/) f[$3]++; \
else if(/Mv-cat-C212-A1:[0-9]+\.[0-9]+,Mv-sup-1065-A1/ || /Mv-cat-C212-A2:[0-9]+\.[0-9]+,Mv-sup-1065-A2/) p[$3]++; \
else if(/Mv-cat-C212-A1:[0-9]+\.[0-9]+,Mv-cat-C212-A2/ && /Mv-sup-1065-A1:[0-9]+\.[0-9]+,Mv-sup-1065-A2/) s[$3]++}END{for(i in f) print "full",i,f[i]; \
for(j in p) print "partial",j,p[j]; for(k in s) print "species",k,s[k]}' g2ogf2stratum2tree_reordered \
| sort -k2,2 -k1,1