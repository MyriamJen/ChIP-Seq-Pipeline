# ChIP-Pipeline

## 1. Checking the sequence quality by using FastQC
All samples are run in a loop for the initial FastQC analysis. For this a file called allsamples.txt is created, where all the files are named 
without their suffix (to circumvent strange names coming out of the loop). Their suffix was added to their variable name and their output directed 
to a new folder "1_FastQC" that branches from the parent directory holding the files. The quality assessment was done manually.

```bash
for i in $(cat allsamples.txt);
do fastqc $i".fq.gz" > 1_FastQC/$i;
done
```

## 2. Trim the sequences by TrimGalore
This step is performed to ensure analysis of high quality reliable reads only. In this step, first Cutadapt is used to trim off any remaining adapter
 sequences and to clear the files from highly overrepresented sequences. Then FastQC is used again to control the results. Note, that this treatment 
was sufficient for all samples except for the ring stage parasite samples sequenced at AGRF in 2018, where several steps were needed to clear off all 
overrepresented sequences (TrimGalore1 (2) --> Cutadapt1 (2a) --> TrimGalore2 (2b) --> Cutadapt2 (2c)). 
This time the allsample.txt file was again generated from the names of all files, just without their number, which indicates their reading direction so
 that the paired reads could be analysed in TrimGalore together. This time the output was directed into the "2_TrimGalore" folder and the files were 
ended with "_val_1" or "_val_2" automatically. Additionally the FastQC results were directly put into "3_Trim_FastQC" folder. Also, 4 cores were used (j), 
all reads that were <20 bases after trimming were removed (length), paired end reads were taken into account (paired) and the illumina adapters were 
searched for (illumina).

```bash
for i in $(cat allsample.txt);
do TrimGalore $i"_1.fq.gz" $i"_2.fq.gz" -o 2_TrimGalore --fastqc_args "--outdir 3_Trim_FastQC" \
-j 4 --illumina --paired --length 20;
done
```

## 3. Aligning the sequences using Bowtie2
Subsequently, the reads were aligned to the reference genome version 28 of Plasmodium falciparum 3D7. The reference genome was obtained from the PlasmoDB website. For this first, the index file was generated. Then again the allname.txt file was used (the one used for TrimGalore) and the resulting files from TrimGalore were used. The output, samfiles and log files, were directed into "4_Samfiles" folder. The mapping quality is stored in column 5 of the resulting samfile and can be plotted as a histogram to get an idea of the distribution and a suitable quality cutoff later on. The values are extracted and sorted using the bash shell first (the '-out.txt' and '-out-sorted.txt' files can be removed afterwards), then the histograms are plotted using R. The resulting files are stored in the "4a_Mapping_Quality" folder. The histograms are saved in a separate R markdown "Mapping_Quality_Graphs_output.Rmd". Here only one file is shown as an example. Overall, the quality seemed to be quite good (usually mostly 42) but there were also reads with quality 0 and especially in the AGRF samples, most samples showed a small peak around 24. Therefore, for the generation of bw files later on, a quality cutoff of 28 was chosen.

```bash
#Generation of the index files
bowtie2-build --threads 8 -f PlasmoDB-28_Pfalciparum3D7_Genome.fasta index
#Aligning the reads to the indexed genome
for i in $(cat allname.txt); do
bowtie2 -x index -1 2_TrimGalore/$i'_R1_001_val_1.fq.gz' \
-2 2_TrimGalore/$i'_R2_001_val_2.fq.gz' \
-S 4_Samfiles/$i'.sam' -t -p 8 2> 4_Samfiles/$i'.log';
done
#other options:
#-5 *number* (trim a number of bases from the 5' end)
#-3 *number* (trim a number of bases from the 3' end)
#--phred33 or phred64 (depending on sequencer)
#-l *number* (minimum distance between paired reads, default 0)
#-X *number* (maximum distance between paired reads, default 500)
#--no-mixed (disables single reads if no alignment pair is found)
#-t (write time to output)
#-p *number* (number of threads to be used)
#Visualize the mapping quality using bash
for i in $(cat allnames.txt); do
awk 'BEGIN {FS=OFS="\t"} {print $5}' $i'.sam' > $i'-out.txt' | sort -n $i'-out.txt' | uniq -c > $i'-hist.txt';
done

#in case you already deleted the sam files:
for i in $(cat allnames.txt); do
samtools view 5_Bamfiles/$i'.bam' | awk 'BEGIN {FS=OFS="\t"} {print $5}' | sort -n | uniq -c > $i'-hist.txt';
done
```
```R
#First read the histogram table and name the variable
GD6_L1_1 <- read.table("/Volumes/NGS_1/ChIP-Seq/Run_2015_Shamista/D6/D6_2/4_Samfiles/4a_Mapping_Quality/GD6-L1_S9_L001-hist.txt", quote="\"", comment.char="")
#Then plot as barplot
barplot(GD6_L1_1$V1, width = 1, space = 0.1, names.arg = GD6_L1_1$V2, main = "GD6_L1_1", xlab = "Mapping Quality", ylab = "Frequency")
```

## 4. Transform the samfiles in bamfiles
In the next step, the samfiles are transformed into their binary version, the bamfiles. These are then indexed. Then also replicates are merged, for some files of Run_2015_Shamista and Run_2018_AGRF, which were run on multiple lanes. For all this, samtools is used. For the loop we again need the allname.txt file and access the samfiles in folder "4_Samfiles". The resulting bamfiles are stored in teh "5_Bamfiles" folder. If files are merged, the resulting files will be stored in the temporary "5a_Bamfiles_merged" folder and subsequently sorted and indexed again and then stored in the final folder "5b_Bamfiles_merged_sorted". For merging, again a file "name.txt" is needed, which states the names of the files to be merged, only the lane information is then added as suffix to specify the files.

```bash
for i in $(cat allname.txt); do
samtools view -Sb 4_Samfiles/$i'.sam' | samtools sort - > 5_Bamfiles/$i'.bam' && \
samtools index 5_Bamfiles/$i'.bam';
done
#Merging Bamfiles
for i in $(cat name.txt); do
samtools merge -f 5a_Bamfiles_merged/$i'.bam' 5_Bamfiles/$i'_L001.bam' 5_Bamfiles/$i'_L002.bam' 5_Bamfiles/$i'_L003.bam';
done
for i in $(cat name.txt); do
samtools sort 5a_Bamfiles_merged/$i'.bam' > 5b_Bamfiles_merged_sorted/$i'_sorted.bam' && \
samtools index 5b_Bamfiles_merged_sorted/$i'_sorted.bam';
done
```

## 5. Generate Bigwigfiles for visualisation
In a further step, the bam files are transformed into bigwigfiles, that can be visualised in IGV or used for further analysis of the coverage data. The package used for this purpose is deeptools and the subcommands are bamCoverage (for "normal" bw files) and bamCompare (for logratio bw files). Here basically, the coverage per bin is calculated and can be normalised in different ways: RPKM (reads per kilobase per million mapped fragments), CPM (counts per million), BPM (bins per million mapped reads), RPCG (1x depth - reads per genome coverage). If the normalisation option is used, scaleFactorsMethod has to be set to 'None', otherwise it can also be 'readCount' or 'SES'. The resulting files are stored in "6_Bigwigfiles". Again, the allname.txt file is used for the generation of the "normal" bigwigfiles. For the generation of the logration files, of sample over input (log2), a file with sample names "samples.txt" and a files with the input names "input.txt" is used. Additionally also logratio Bigwigfiles for the merged Bamfiles later on in the Merged Analysis and also Bedgraph files are generated. The bigwigfiles will be needed for the merged analysis with compute_matrix (Part 10.) and the bedgraph files for plotting the UMAP (Part 12.). For the bedgraph files, bs is changed to 300 and the ending is turned into bg instead of bw. The bedgraph files are stored in Merged_Analysis/10_BedGraphfiles_merged_sorted_replicates.

```bash
for i in $(cat allname.txt); do
bamCoverage -p 8 --normalizeUsing 'RPKM' --minMappingQuality 28 --ignoreDuplicates -b 5_Bamfiles/$i'.bam' -o 6_Bigwigfiles/$i'_rpkm_ndq.bw';
done
for i in $(cat samples.txt); do
read j
bamCompare -p 8 --scaleFactorsMethod 'None' --normalizeUsing 'RPKM' --minMappingQuality 28 --ignoreDuplicates -b1 5_Bamfiles/$i'.bam' -b2 5_Bamfiles/$j'.bam' -o 6_Bigwigfiles/$i'-logratio_rpkm_ndq.bw';
done < input.txt
#other options:
#-bs *number* (to set the binsize, default is 50bp)
#--operation *log2 / ratio / substract / add / mean / reciprocal_ratio / first / second* (which operation is used for the comparison of sample and input, default is log2)
#--skipZeroOverZero (skips bins in which input and sample have no coverage)
```

## 6. Generate a MultiBamSummary
For quality control of the samples, a multibamsummary can be generated. This shows if replicates cluster together and if different samples cluster away from each other. This was done for all merged files and finally also for the runs alltogether to see if the data is reporoducable. The input can be as many bam files as wished to be compared, the output is a npz array, that can be used to plot a heatmap or PCA plot. The output files are all stored in the "7_MulitBamSummary" folder. 

```bash
multiBamSummary bins -p 4 --minFragmentLength 50 --maxFragmentLength 3000 -b 5_Bamfiles_merged_sorted/*.bam \
-out 7_MultiBamSummary/*filename*-summary-merge.npz
plotCorrelation -in 7_MultiBamSummary/*filename*-summary-merge.npz -o 7_MultiBamSummary/*filename*-summary-merge-heatmap.pdf \
-c spearman -p heatmap --plotNumbers --plotHeight 42 --plotWidth 29.7 --smartlabels \
--outFileCorMatrix 7_MultiBamSummary/*filename*-summary-merge-heatmap.txt
plotPCA -in 7_MultiBamSummary/*filename*-summary-merge.npz -o 7_MultiBamSummary/*filename*-summary-merge-pca.pdf
```

## 7. Call the peaks
Next, the peak enrichment regions over Input are to be called using MACS2. The results are narrow.Peak files, that can be loaded onto the bw file tracks in IGV to control where the peaks were called. Additionally, the summits are defined in summits.bed files and the peak regions are shown in peaks.xls files. All those resulting files are stored in the "8_Peaks" folder.

```bash
macs2 callpeak -t 5_Bamfiles/*sample*.bam \
-c 5_Bamfiles/*input*.bam --format BAMPE \
 --name *sample*_peaks --outdir 8_Peaks/ --gsize 2.3e7 --nomodel -p 0.001 --call-summits
 
 #other options:
 #--format *BAM / BAMPE* (BAMPE for paired end bam files, otherwise BAM)
 #--gsize 2.3e7 (state the mappable size of the genome)
 #--nomodel (use to have more comparable peak calling between samples by using the same default values instead of individually calculated ones)
 #--mfold *number* *number* (fold change cutoff for peaks, 5 50 default but can be used as 2 50 or 2 100)
 #--extsize *number* (default fragment size d is 200 for nomodel, but can be changed to 147)
 #--min-length *number* (minimum peak length, calculated as a fragment of d)
 #--nolambda (no local correction, only default globally)
 #--qvalue *number* (set q value, default 0.05)
 #--pvalue *number* (set p value, overwrites the q value)
```
