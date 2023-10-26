# Generation of gene expression groups

## 1. Use Multicov to predict heterochromatic areas
In this step, we want to see whether a gene should be considered heterochromatic or not. Per definition, a heterochromatic gene will be covered in H3K9me3.
Therefore, the 500bp before and after the ATG of each gene are extracted by using bedtools multicov. Then the H3K9me3 coverage is analysed in a bivariat 
gaussian model and the probability of the gene to be either heterochromatic or euchromatic is predicted. This is then used to create lists of 
heterochromatic genes and together with expression data is used to sort the genes in groups per stage of

•	heterochromatic
•	silent (FPKM = 0)
•	lowly transcribed (33% lowest FPKMs)
•	moderately transcribed (33% middle FPKMs)
•	highly transcribed (33% highest FPKMs)

Here again for the multicov calculations, the samples and input files are needed. The coverage will be returned from the bam files, the resulting value 
needs to be normalised to the input, therefore it is useful to have both informations in adjacent columns of one file. As input bed file with the desired
regions, the 3D7-genes-atg500-ok.bed file is used, the output files are stored in the “9_Multicov” folder. Finally, the txt files are used in R with the
mixtools package for NormalmixEM. Finally, this information will be used to distribute all genes into the different expression groups - the expression 
data is generated from the RNA-Seq data.

```bash
#Generation of the 3D7-genes-atg500-ok.bed file
##Generate a file with all genes
awk 'BEGIN {FS=OFS="\t"} {if ($3 == "gene") {print $0} else {}}' PlasmoDB-28_Pfalciparum3D7.gff > 3D7-genes.gff
##Loose the apicoplast and mitochondrium
awk 'BEGIN {FS=OFS="\t"} {if ($1 ~ /Pf3D7/) {print $0} else {}}' 3D7-genes.gff > 3D7-genes-noma.gff
##Get the region +/- 500bp around the ATG and output a bed file
awk 'BEGIN {FS=OFS="\t"} {if ($7 == "+") {atg = $4} else {atg = $5} print $1, atg-500, atg+500, ".", ".", $7, $9}' 3D7-genes-noma.gff > 3D7-genes-atg500.bed
##Turn negative values into 1 if present (if gene would start less than 500 away of the chromosome end)
awk 'BEGIN {FS=OFS="\t"} {if ($2 < 1) {print $1, 1, $3, $4, $5, $6, $7} else {print $0}}' 3D7-genes-atg500.bed > 3D7-genes-atg500-ok.bed

#Bedtools Multicov
for i in $(cat samples.txt); do
read j
bedtools multicov -bams 5_Bamfiles/$i'.bam' 5_Bamfiles/$j'.bam' \
-bed 3D7-genes-atg500-ok.bed > 9_Multicov/$i'-atg500.txt';
done < input.txt
```

Use R for NormalmixEM. Here only the code for the day 6 sample from 2015 is shown as an example - this was done for all samples for H3K9me3 and the 
average probability between the replicates was then used further.

```R
library(mixtools)
## mixtools package, version 1.2.0, Released 2020-02-05
## This package is based upon work supported by the National Science Foundation under Grant No. SES-0518772.
#Read the coverage data of H3K9me3 and input ChIP-Seq samples into R

Run_2015_D6_H3K9me3table <- read.delim("/Volumes/NGS_1/Chromatin/ChIP-Seq_Chromatin/Run_2015_Shamista/D6/D6_2/9_Multicov/GD6-L5_S34-atg500.txt", header=FALSE)

#Create column 10 as the ratio between column 8 and 9 (make sure column 9 is not 0)

Run_2015_D6_H3K9me3 <- c(Run_2015_D6_H3K9me3table$V8/Run_2015_D6_H3K9me3table$V9)

#bind the new column 10 to the old table 

Run_2015_D6_H3K9me3tablen <- cbind(Run_2015_D6_H3K9me3table, Run_2015_D6_H3K9me3)

#Define the formula to plot the graphs
plot.normal.components <- function(mixture,component.number) { curve(mixture$lambda[component.number] * 
                                                                       dnorm(x,mean=mixture$mu[component.number], 
                                                                             sd=mixture$sigma[component.number]), add=TRUE)
}

#Use NormalmixEM to find the probabilities for euchromatic vs heterochromatic genes

OutputRun_2015_D6_H3K9me3 <- normalmixEM(Run_2015_D6_H3K9me3tablen$Run_2015_D6_H3K9me3, lambda = NULL, mu = NULL,
                                         sigma = NULL, k = 2, mean.constr = NULL, sd.constr = NULL, epsilon = 1e-05,
                                         maxit = 1000, maxrestarts = 1000, verb = TRUE, fast = TRUE, ECM = TRUE,
                                         arbmean = TRUE, arbvar = TRUE)
                                         
summary(OutputRun_2015_D6_H3K9me3)

#Plot the distribution and the fitted curves

plot(hist(Run_2015_D6_H3K9me3tablen$Run_2015_D6_H3K9me3,breaks=100),col="green",border="black",freq=FALSE, prob=TRUE,
     xlab="Readtags over Input",main="Histogram of H3K9me3 Readtags")
lines(density(Run_2015_D6_H3K9me3tablen$Run_2015_D6_H3K9me3),lty=2) 

sapply(1:2,plot.normal.components,mixture=OutputRun_2015_D6_H3K9me3)
 
#Get the results from NormalmixEM into the data table

Run_2015_D6_H3K9me3_HC <- cbind(Run_2015_D6_H3K9me3tablen, OutputRun_2015_D6_H3K9me3$posterior)

#Save the table with the HC information
write.table(Run_2015_D6_H3K9me3_HC, "Run_2015_D6_H3K9me3_HC.txt", sep="\t", row.names = F, quote = F, col.names = F)
```

## 2. Merge the HC information with Expression Data
For this part, the RNA-Seq data is needed. To determine the groups described, a file with the HC probability and the expression data for every stage per 
gene is needed. When the right HC, expression, and the gene file are chosen, they can be fused with filemerge.awk. The resulting groups are stored in the 
folder “7_Genesubsets” in the RNA-Seq section.

```awk
#!/bin/awk -f

# the number of phases depends on the number of files red, first the transcript.gtf files are red, then the HC files, then the 3D7 reference file. 

BEGIN {
    FS=OFS="\t"
    phase = 1
}
ENDFILE {
    phase += 1
}

function get_info_value(info, key, sep, isQuoted) {
    split(info, arr, ";")
    # Entries are formattel like:
    # name blank value
    # where the value is enclosed in double quotes
    if (sep == "") {
        sep = " "
    }
    if (isQuoted == "") {
        isQuoted = 1
    }
    for (entry in arr) {
        split(arr[entry], e, sep)
        if (e[1] == key) {
            if (isQuoted == 1)
                return substr(e[2], 2, length(e[2]) - 2)
            else
                return e[2]
        }
    }
    return "NOT FOUND"
}

phase < 5 {
    info = $9
    id = get_info_value(info, "gene_id")
    FPKM = get_info_value(info, "FPKM")
    genarr[id,phase] = FPKM
}

phase >= 5 && phase <8 {
    info = $9
    id = get_info_value(info, "ID", "=", 0)
    HC = $12
    genarr[id,phase] = HC
}

phase == 8 && $3 == "gene" {
    info = $9
    id = get_info_value(info, "ID", "=", 0)
    if (genarr[id,1] != "") {
        FPKM1 = genarr[id,1]
        FPKM2 = genarr[id,2]
        FPKM3 = genarr[id,3]
        FPKM4 = genarr[id,4]
        HC1 = genarr[id,5]
        HC2 = genarr[id,6]
        HC3 = gerarr[id,7]
        print $1, $4, $5, ".", ".", $7, id, FPKM1, FPKM2, FPKM3, FPKM4, (FPKM1 + FPKM2 + FPKM3 + FPKM4)/4, HC1, HC2, HC3, (HC1 + HC2 + HC3)/3
    }
}
```

The resulting file will be in bed format with extra columns for FPKM, WSKH and WSKE.

```bash
awk -f filemerge.awk transcripts.gtf *_HC.txt 3D7-genes-noma.bed > *_match.bed
```

Then divide the genes in the desired groups:

```bash
#!/bin/bash
#here it is assumed that the average FPKM is in column 8, and WSKH in column 9)
#Get files with only heterochromatic genes
for i in $(cat names.txt); do
    awk 'BEGIN {FS=OFS="\t"} {if ($9 > 0.99999) {print $0} else {}}' $i'.bed' > $i'_HC.bed' && echo $i "is ready";
done

#Get files with only silent but not heterochromatic genes
for i in $(cat names.txt); do
    awk 'BEGIN {FS=OFS="\t"} {if ($8 == 0 && $9 <= 0.99999) {print $0} else {}}' $i'.bed' > $i'_silent.bed' && echo $i "is ready";
done

#Get files with only transcribed genes
for i in $(cat names.txt); do
    awk 'BEGIN {FS=OFS="\t"} {if ($8 != 0 && $9 <= 0.99999) {print $0} else {}}' $i'.bed' > $i'_tr.bed' && echo $i "is ready";
done

#Sort the files containing the transcribed genes in -r reverse order (highest value on top), -n numerical, -k 8 sort for values in column 8 (thus RPKM values)
for i in $(cat names.txt); do
    sort -r -n -k 8 $i'_tr.bed' > $i'_tr_sorted.bed' && echo $i "is ready";
done

#Create the transcription level subgroup files
for i in $(cat names.txt); do
    cat $i'_tr_sorted.bed' | wc -l > $i'_num.txt' &&
    awk '{print int($1/3)}' $i'_num.txt' > $i'_num3.txt' &&
    awk '{print int($1*2)}' $i'_num3.txt' > $i'_num32.txt' &&
    echo $i "is ready";
done

for i in $(cat names.txt); do
    head -n $(cat $i'_num3.txt') $i'_tr_sorted.bed' > $i'_tr_top.bed' && echo $i "is ready" && \
    tail -n $(cat $i'_num3.txt') $i'_tr_sorted.bed' > $i'_tr_bottom.bed' && echo $i "is ready" && \
    head -n $(cat $i'_num32.txt') $i'_tr_sorted.bed' | tail -n $(cat $i'_num3.txt') > $i'_tr_middle.bed' && echo $i "is ready";
done
