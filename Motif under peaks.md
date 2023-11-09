# Analysis of DNA motifs at ChIP-Seq peaks

Column 2 and 3 are the coordinates of the summit, column 5 is the fold change in the summits.bed files

```{bash}
#--> First exclude telomeres (use before first gene and after last gene as borders), then use the top 500 summits

awk 'BEGIN {FS=OFS="\t"} {if ($1 == "Pf3D7_01_v3" && $2 > 30000 && $3 < 615000) {print $0} else \
{ if ($1 == "Pf3D7_02_v3" && $2 > 25000 && $3 < 925000) {print $0} else \
{ if ($1 == "Pf3D7_03_v3" && $2 > 37000 && $3 < 1040000) {print $0} else \
{ if ($1 == "Pf3D7_04_v3" && $2 > 25000 && $3 < 1180000) {print $0} else \
{ if ($1 == "Pf3D7_05_v3" && $2 > 20000) {print $0} else \
{ if ($1 == "Pf3D7_06_v3" && $2 > 600 && $3 < 1385000) {print $0} else \
{ if ($1 == "Pf3D7_07_v3" && $2 > 20000 && $3 < 1425000) {print $0} else \
{ if ($1 == "Pf3D7_08_v3" && $2 > 20000 && $3 < 1445000) {print $0} else \
{ if ($1 == "Pf3D7_09_v3" && $2 > 20000 && $3 < 1505000) {print $0} else \
{ if ($1 == "Pf3D7_10_v3" && $2 > 30000 && $3 < 1650000) {print $0} else \
{ if ($1 == "Pf3D7_11_v3" && $2 > 25000 && $3 < 2035000) {print $0} else \
{ if ($1 == "Pf3D7_12_v3" && $2 > 20000 && $3 < 2250000) {print $0} else \
{ if ($1 == "Pf3D7_13_v3" && $2 > 20000 && $3 < 2895000) {print $0} else \
{ if ($1 == "Pf3D7_14_v3" && $2 > 1500 && $3 < 3258000) {print $0} else {}}}}}}}}}}}}}}}' summits.bed > summits_notel.bed #example file!

#sort by FC and get 500 highest
sort -r -nk5 summits_notel.bed | head -n 500 > top500_summits.bed #example file!

#get +/- 50bp around summit
awk 'BEGIN {FS=OFS="\t"} {print $1, $2-49, $3+49}' top500_summits.bed > top500_summits_100bp.bed #example file!

#get the sequence
bedtools getfasta -fi PlasmoDB-59_Pfalciparum3D7_Genome.fasta -bed top500_summits_100bp.bed -fo top500_summits_100bp.fa #example file!

#find the motif
meme top500_summits_100bp.fa -dna -nmotifs 12 -evt 0.05 -p 8 -V -minw 4 -maxw 12 -oc SAMPLENAME -dna -revcomp; done
```
