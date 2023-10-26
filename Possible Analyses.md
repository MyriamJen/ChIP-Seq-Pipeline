# 1. Correlation between samples using multiBigwigSummary

```{bash}
multiBigwigSummary bins -p 4 -b 'bigwigfiles' -bs 150 -v -out result.npz

plotCorrelation -in results.npz -o results.svg -c spearman -p scatterplot --skipZeros -l 'label' --plotHeight 4 --plotWidth 5 
```
## Comments:
-p number of processors
-b insert bigwig files here, two or more are possible
-bs bin size to use for calculations, must be more or equal to bin size used for generation of bw files
-c choose between pearson and spearman correlation
-p plottype
-l insert one label per bw file


# 2. Create Heatmaps and Lineplots using deeptools computeMatrix

```{bash}
computeMatrix scale-regions -p "max" \
-R bedfiles \
-S bigwigfiles \
-o results \
-m 2000 -b 2000 -a 2000 --startLabel "ATG" --endLabel "STOP" --sortRegions keep --verbose \
--samplesLabel "label"

plotHeatmap -m GD6F_all_new -o GD6F_all_heatmap_new.svg \
--plotType lines \
--sortRegions keep \
--colorMap seismic \
--whatToShow "heatmap and colorbar" \
--startLabel "ATG" --endLabel "STOP" --regionsLabel "label" \
--samplesLabel "label" --zMin -2 --zMax 2
```

## Comments:
-m size of region in the plot
-a number of bases ups of region
-b number of bases downstream of region
--startLabel labels the start of the plotted region
--endLabel labels the end of the plotted region
--sortRegions you can choose from keep, ascending, descending
--samplesLabel label all the bigwigfiles, put labels in quotation marks

--plotType choose what type of lines (lines, sd etc)
--sortRegions you can choose from keep, ascending, descending
--colorMap choose what type of colour scheme you want, seismic is red/blue
--whatToShow choose from heatmap, lineplot, colorbar or combinations of the three
--regionsLabel label all the bedfiles for the regions
--zMin/zMax set the minimum and maximum of the scalebar for z (in heatmap), same also possible for y (in lineplot)

