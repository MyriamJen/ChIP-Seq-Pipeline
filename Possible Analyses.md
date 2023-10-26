1. Correlation between samples using multiBigwigSummary

multiBigwigSummary bins -p 4 -b 'bigwigfiles' -bs 150 -v -out result.npz

plotCorrelation -in results.npz -o results.svg -c spearman -p scatterplot --skipZeros -l 'label' --plotHeight 4 --plotWidth 5 

Comments:
-p number of processors
-b insert bigwig files here, two or more are possible
-bs bin size to use for calculations, must be more or equal to bin size used for generation of bw files
-c choose between pearson and spearman correlation
-p plottype
-l insert one label per bw file


2. Create Heatmaps and Lineplots using deeptools computeMatrix


