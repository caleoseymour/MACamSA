library('Rmisc')
library('phyloseq')
library('dplyr')
library('Rmisc')
library('ggplot2')
library("grid")
library("RColorBrewer")
library("vegan")

## Function used to generate the seed at the time of execution.
# set.seed(as.integer(runif(1)*2e9))
# rseed = .Random.seed
# seed_date = date()
# save(".Random.seed",file="random_seed.RData")

## For reproducibility
load("random_seed.RData")
rseed = .Random.seed
set.seed(rseed)

source('importTSVPhyseq.R')
physeq = importPhyseqFromQiime2(otutable = '../exported/plaintext_feature_table.tsv',
                                samdata = '../metadata.csv',
                                taxtable = '../exported/silva-taxonomy.tsv')
source('makeBarGraphs.R')
f_physeq = prune_samples(rownames(sample_data(physeq))[sample_data(physeq)$trmnt_type_s == "CamSa"],physeq)
f_physeq = prune_taxa(taxa_sums(f_physeq) > 0, f_physeq)
## Make some fixes to the taxonomy...
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Eubacterium] coprostanoligenes group',5] = 'D_4__Eubacteriaceae'
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Eubacterium] coprostanoligenes group',6] = 'D_5__Eubacterium'
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Eubacterium] xylanophilum group',5] = 'D_4__Eubacteriaceae'
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Eubacterium] xylanophilum group',6] = 'D_5__Eubacterium'
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Eubacterium] ruminantium group',5] = 'D_4__Eubacteriaceae'
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Eubacterium] ruminantium group',6] = 'D_5__Eubacterium'
tax_table(f_physeq)[tax_table(f_physeq)[,6] == 'D_5__[Ruminococcus] torques group',6] = 'D_5__Blautia'

plot_bar_ext(f_physeq,fill='Genus',cutoff=0.05)
plot_bar_ext(f_physeq,fill='Family',cutoff=0.05)
plot_bar_ext(f_physeq,fill='Order',cutoff=0.05)
plot_bar_ext(f_physeq,fill='Class',cutoff=0.05)
plot_bar_ext(f_physeq,fill='Phylum',cutoff=0.05)

#r_physeq = rarefy_even_depth(physeq)

#source('alphaDiversityPlots.R')
source('adonisAndSimper.R')
source('Lactobacillus-Bifidobacterium_Search.R')
