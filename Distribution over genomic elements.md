# Distribution over genomic elements

To learn more about the distribution of the called ChIP-Seq peaks, you can look at their distribution over the genome. You can distinguish different genomic elements: coding regions (exons+introns), promoter containing regions (divergent and tandem regions) and regions downstream of genes (convergent regions).

For this, you need to create bed files with the information on the coordinates of the regions. Use the version of the genome that you used for aligning the data.

## 1. Create files with exon and intron information

```{R}
library(GenomicFeatures)
library(rtracklayer)

gtf <- makeTxDbFromGFF("PATH/PlasmoDB-59_Pfalciparum3D7.gff")
exons <- exonsBy(gtf, by="tx")

exons <- reduce(exons)
exons <- exons[sapply(exons, length) > 1]

export(exons, "PATH/exons.gff", format = "gff3")

#make introns
introns <- lapply(exons, function(x) {
    gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
        end=max(end(x))), 
        strand=strand(x)[1])
    db = disjoin(c(x, gr))
    ints = db[countOverlaps(db, x) == 0]
    if(as.character(strand(ints)[1]) == "-") {
        ints$exon_id = c(length(ints):1)
    } else {
        ints$exon_id = c(1:length(ints))
    }
    ints
})
intron <- introns[sapply(introns, length) > 1]
intronn <- GRangesList(intron)

export(intronn, "PATH/introns.gff", format = "gff3")
```
```{bash}
#use awk to make bedfiles
awk 'BEGIN {FS=OFS="\t"} {if ($3 == "exon") {print $1, $4, $5, $6, $8, $7, $9} else {}}' introns.gff > introns.bed
awk 'BEGIN {FS=OFS="\t"} {if ($3 == "exon") {print $1, $4, $5, $6, $8, $7, $9} else {}}' exons.gff > exons.bed
```

## 2. Create intergenic regions and then sort them
```{bash}
#generate genes
awk 'BEGIN {FS=OFS="\t"} {if ($3 ~ /gene/) {print $1, $4, $5, $6, $8, $7, $9} else {}}' PlasmoDB-59_Pfalciparum3D7.gff > PlasmoDB-59_Pfalciparum3D7_genes.bed
bedtools sort -i PlasmoDB-59_Pfalciparum3D7_genes.bed > PlasmoDB-59_Pfalciparum3D7_genes_sorted.bed

#create the intergenic regions (per chromosome) --> do for all chromosomes
awk 'BEGIN {FS=OFS="\t"} {if ($1 == "Pf3D7_01_v3") {print $0} else {}}' PlasmoDB-59_Pfalciparum3D7_genes_sorted.bed > PlasmoDB-59_Pfalciparum3D7_genes_sorted_chr1.bed

#start of the intergenic regions is column 2 minus 1 (shifted one line down), direction of the previous gene is column 6 shifted one line down, end is column 3 plus 1, direction of the next gene is column 6 --> do that for all 14 chromosomes
awk 'BEGIN {FS=OFS="\t"} NR>1 {print last} {last=$0}' PlasmoDB-59_Pfalciparum3D7_genes_sorted_chr1.bed | awk '{print $1, $2-1, $6}' > PlasmoDB-59_Pfalciparum3D7_intergenes_start_sorted_chr1.bed #cut off last line to leave out the right telomere
awk 'BEGIN {FS=OFS="\t"} NR>1 {print $1, $3+1, $6}' PlasmoDB-59_Pfalciparum3D7_genes_sorted_chr1.bed > PlasmoDB-59_Pfalciparum3D7_intergenes_end_sorted_chr1.bed #start at second line to leave out the left telomere
paste PlasmoDB-59_Pfalciparum3D7_intergenes_start_sorted_chr1.bed PlasmoDB-59_Pfalciparum3D7_intergenes_end_sorted_chr1.bed > PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr1.txt
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $5, ".", ".", $3, $6}' PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr1.txt > PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr1.bed

cat PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr1.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr2.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr3.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr4.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr5.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr6.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr7.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr8.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr9.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr10.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr11.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr12.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr13.bed PlasmoDB-59_Pfalciparum3D7_intergenes_sorted_chr14.bed > PlasmoDB-59_Pfalciparum3D7_intergenes_sorted.bed 

#define the subgroups: tandem (++ or --), divergent (-+), convergent (+-)
awk 'BEGIN {FS=OFS="\t"}{if ($6 == "+" && $7 == "+") {print $0} else {if ($6 == "-" && $7 == "-") {print $0} else {}}}' PlasmoDB-59_Pfalciparum3D7_intergenes_sorted.bed > PlasmoDB-59_Pfalciparum3D7_tandem.bed

awk 'BEGIN {OFS="\t"}{if ($6 == "-" && $7 == "+") {print $0} else {}}' PlasmoDB-59_Pfalciparum3D7_intergenes_sorted.bed > PlasmoDB-59_Pfalciparum3D7_divergent.bed

awk 'BEGIN {OFS="\t"}{if ($6 == "+" && $7 == "-") {print $0} else {}}' PlasmoDB-59_Pfalciparum3D7_intergenes_sorted.bed > PlasmoDB-59_Pfalciparum3D7_convergent.bed
```

## 3. Filter your peaks 
so that you have only one peak per peak regions instead of multiple due to different summits

```{bash}
#work with unique peaks instead of one peak with many summits
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, ".", ".", "."}' PfBDP4loxPHA_R_1_peaks_peaks.narrowPeak | sort -u > PfBDP4loxPHA_R_1.bed #example file!
```

## 4. Intersect the peaks with the genetic elements to see their distribution
Note: the sum of the peaks intersecting the elements is equal or higher than the number of peaks, as peaks can intersect multiple regions.

```{bash}
#intersect with exons (use uniq as it is only interesting if the peak has any hit or not)
bedtools intersect -u -a PfBDP4loxPHA_R_1.bed -b exons.bed | wc -l  #example file!

#intersect with introns
bedtools intersect -u -a PfBDP4loxPHA_R_1.bed -b introns.bed | wc -l  #example file!

#intersect with convergent regions
bedtools intersect -u -a PfBDP4loxPHA_R_1.bed -b PlasmoDB-59_Pfalciparum3D7_convergent.bed | wc -l  #example file!

#intersect with divergent regions
bedtools intersect -u -a PfBDP4loxPHA_R_1.bed -b PlasmoDB-59_Pfalciparum3D7_divergent.bed | wc -l  #example file!

#intersect with tandem regions
bedtools intersect -u -a PfBDP4loxPHA_R_1.bed -b PlasmoDB-59_Pfalciparum3D7_tandem.bed | wc -l  #example file!
```






