# Detecting islands of differentiation in nightingale genome

The goal of the analysis is to find regions the genomes of Common and Thrush
nightingales, which could contribute to the recent speciation event. We used
SNP variants, per site Fst and Dxy as indicators of possible speciation,
alignment to Zebra finch as  a proxy for genomic location of the variants,
floating windows to scan along chromosomes, and bootstrap to detect
significanly higher values.

## Variant filtering 
To get more reliable results, variants should be filtered on some biological
significance. Filter out nightingale contigs mapped to too wide range in the
zebra finch, filter out variants on some quality threshold, for Dxy filter
variants with too few called  individuals.

## Window calculation
Cover all chromosomes with equally spaced intervals. Use at least 1000 intervals
for the longest chromosome (?). Most of the calculations is done in pure `data.frame`s,
with the help of `dplyr`, `tidyr` and visualised with `ggplot2`.
Now it's using `GenomicRanges` package only to find the variants belonging to each window. 
I'd like to get rid of any BioConductor traces from the whole project, substituting
it with `data.table::foverlaps`.

## Bootstrap
Bootstrap was performed via permuting the values (*Fst*, *Dxy*) assigned to each genetic 
variant within chromosomes, the variants remain at original positions. We sampled
25,000 permutations, taking the maximum of the windowed measure at given grid point
as the threshold value for 'maximal Fst/Dxy' arising by pure chance in the given variants 
and chromosomal arrangement.

## Islands of speciation
Areas where both Fst and Dxy exceeded the threshold found in the bootstrap
were taken  as putative 'islands of speciation'. We used `biomaRt` to get
information on genes found in those  regions.

## Gene set enrichment
After trying several online tools and having non-reproducibility issues, we decided to 
run the analysis on simple KEGG pathways on our own. Using simple gene lists in `data.frame`s
and Fisher's exact and Binomial tests.

## TODO
- more details on the analysis

## old notes from the US
Dxy windows
dxy - scaled by contig size

windowed fixed diffs / (fixed + shared) should be higher where fst and dxy are both high

dnds over several birds
NUP genes in drosophila (daven presgraves)

paper:
- one or two? (Z and selection over genome/lineages)
- put all informative figures together, decide how many papers;)



