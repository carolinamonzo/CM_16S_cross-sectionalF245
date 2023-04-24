# load libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(phangorn); packageVersion("phangorn")
library(DECIPHER); packageVersion("DECIPHER")
library(vegan); packageVersion("vegan")
library(DESeq2); packageVersion("DESeq2")
library(tidyr); packageVersion("tidyr")
library(svglite)
library(pairwiseAdonis)
library(biomformat)
library(dplyr)
library(reshape2)
library(argparser)
#library(microbiomeSeq)
theme_set(theme_bw())

path <- "~/workspace/16S_final/CM_16S_cross-sectionalF245/"

new_day <- gsub("-", "", as.character(Sys.Date()))

## Make parser
parser <- arg_parser("Script to generate alpha and beta diversity values and plots")

# Add command line arguments
parser <- add_argument(parser, "--path", help = "Path to project")
parser <- add_argument(parser, "--nochim", help = "Path to SeqtabNochim from error_and_cleanTables.R")
parser <- add_argument(parser, "--taxa", help = "Path to Taxa from error_and_cleanTables.R")
parser <- add_argument(parser, "--metadata", help = "Path to metadata file")
parser <- add_argument(parser, "--count_tab", help = "Path to counts_merged from error_and_cleanTables.R")

args <- parse_args(parser)

path <- args$path

#################
### Read our files for analysis
#################
seqtab.nochim <- readRDS(args$nochim)
taxa <- readRDS(args$taxa)
sam <- read.table(args$metadata, header = T, row.names = 1, check.names = F, sep = ";")
count_tab <- read.table(args$count_tab, header = T, row.names = 1, check.names = F, sep = "\t")

#################
### Use our normal files
#################

seqtab.nochim <- readRDS(paste0(path, "analysis/seqtab_merge3/CLEAN_merged_seqtabNochim_20210215.rds"))
taxa <- readRDS(paste0(path, "analysis/seqtab_merge3/CLEAN_taxonomy_merged_20210615.rds"))

sam <- read.table(paste0(path, "metadata/metadata_F245_16S.csv"), header=T, row.names=1, check.names=F, sep=";")
count_tab <- read.table(paste0(path, "analysis/seqtab_merge3/mergedQC/CLEAN_ASVs_counts_merged_20210615.tsv"), 
                        header = T, row.names = 1, check.names = F, sep = "\t")
#################
### Make directories
#################
dir.create(paste0(path, "analysis/intermediate"))
dir.create(paste0(path, "analysis/plots"))
dir.create(paste0(path, "analysis/plots/alpha_beta_diversity"))
#################
### Make tree and put into phyloseq
#################
# Make the phylogenetic tree
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences
# Run sequence alighment using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)
# Change sequence alignment outpu to phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# Make distance matrix
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
# Internal maximum likelihood
fit = pml(treeNJ, data = phang.align)
# Change negative edges length to 0
fitGTR <- update(fit, k=4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
sample_data(sam),
tax_table(taxa), phy_tree(fitGTR$tree))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
ps <- prune_samples(sample_sums(ps) >= 1, ps)
taxa_names_ps <- paste0("ASV", seq(ntaxa(ps)))

# Root tree
set.seed(123)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root =TRUE)
# Check that its rooted
is.rooted(phy_tree(ps))

saveRDS(ps, paste0(path, "analysis/seqtab_merge3/phyloseq_obj_", new_day, ".rds"))

ps <- readRDS(paste0(path, "analysis/seqtab_merge3/phyloseq_obj_20210625.rds"))

#################
### Alpha diversity
#################
write.table(estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson")), paste0(path, "analysis/intermediate/alpha_values_", new_day, ".csv"), sep = ";", quote = F)

pdf(paste0(path,"analysis/plots/alpha_beta_diversity/alpha_Shannon_RGrid_", new_day, ".pdf"))
p = plot_richness(ps, color = "Treatment", measures = c("Shannon"), x = "Months") + theme(legend.title = element_blank()) + geom_point()
p +  scale_color_manual(values = c("magenta", "darkgreen", "gold", "dodgerblue", "red")) + facet_wrap(~Treatment)
dev.off()

pdf(paste0(path,"analysis/plots/alpha_beta_diversity/alpha_Simpson_RGrid_", new_day,".pdf"))
p = plot_richness(ps, color = "Treatment", measures = c("Simpson"), x = "Months") + theme(legend.title = element_blank()) + geom_point()
p +  scale_color_manual(values = c("magenta", "darkgreen", "gold", "dodgerblue", "red")) + facet_wrap(~Treatment)
dev.off()

#################
### UniFraq
#################
uf_distance  <- phyloseq::distance(ps, "uUniFrac")
wuf_distance <- phyloseq::distance(ps, "wUniFrac")
# Save distance matrices for python
write.table(as.data.frame(as.matrix(uf_distance)), paste0(path, "analysis/intermediate/unifraq_unweighted_dist_", new_day, ".csv"), sep = ";", quote = F)
write.table(as.data.frame(as.matrix(wuf_distance)), paste0(path, "analysis/intermediate/unifraq_weighted_dist_", new_day, ".csv"), sep = ";", quote = F)
# Add meta for analysis
df_d <- melt(as.matrix(uf_distance), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
write.table(df_d_new, file = paste0(path, "analysis/intermediate/unifraq_unweighted-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
df_d <- melt(as.matrix(wuf_distance), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
write.table(df_d_new, file = paste0(path, "analysis/intermediate/unifraq_weighted-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
# Plot
uf_ord <- ordinate(ps, method="PCoA", distance=uf_distance)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/unifraq_unweighted_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, uf_ord, color="Treatment") + theme(aspect.ratio=1) +  scale_color_manual(values = c("magenta", "darkgreen", "gold", "dodgerblue", "red"))
dev.off()

wuf_ord <- ordinate(ps, method="PCoA", distance = wuf_distance)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/unifraq_weighted_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, wuf_ord, color="Treatment") + theme(aspect.ratio=1) +  scale_color_manual(values = c("magenta", "darkgreen", "gold", "dodgerblue", "red"))
dev.off()

### Extra for Friday seminar presentation ###
comp = subset_samples(ps, Treatment != "AL_DR12M")
#comp = subset_samples(comp, Treatment != "AL_DR16M")
#comp = subset_samples(comp, Treatment != "AL_DR20M")
comp = subset_samples(comp, Treatment != "AL_lifelong")
comp = subset_samples(comp, Treatment != "DR_lifelong")
comp = subset_samples(comp, Months == 24)
wuf_distance_small <- phyloseq::distance(comp, "uUniFrac")
wuf_ord <- ordinate(ps, method="PCoA", distance = wuf_distance_small)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/5M_ALDRSbray_PCoA_TreatColor_", new_day, ".pdf"), width = 4, height = 2.5)
plot_ordination(ps, wuf_ord, color="Treatment") + theme(aspect.ratio=1, axis.text=element_text(size=32), axis.title=element_text(size=32), 
              plot.title = element_text(size=32)) + geom_point(size = 2.5) + scale_alpha_manual(values = c(0, 0, 1, 1)) +  scale_color_manual(
              values = c("dodgerblue", "red")) + stat_ellipse() + theme_classic() 
              
dev.off()

# Pairwise adonis of this:
print(pairwise.adonis(wuf_distance_small, sample_data(comp)$Treatment, p.adjust.m = "BH"))

#################
### Beta bray curtis
#################
bray_dist = phyloseq::distance(ps, method="bray", weighted=F)

# Store beta diversity for analysis
df_d <- melt(as.matrix(bray_dist), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
# Save the dataframe to plot in python
write.table(df_d_new, file = paste0(path, "analysis/intermediate/beta_values-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
# Save the matrix
write.table(as.data.frame(as.matrix(bray_dist)), paste0(path, "analysis/intermediate/beta_bray_values_", new_day, ".csv"), sep = ";", quote = F)

# Continue
bray_ord = ordinate(ps, method="PCoA", distance=bray_dist)
#plot_ordination(ps, bray_ord, color="Treatment", shape = "Age_char") + theme(aspect.ratio=1) +  scale_color_manual(values = c("darkgreen", "gold", "dodgerblue", "red", "black"))
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/beta_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, bray_ord, color="Treatment") + theme(aspect.ratio=1, axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=18)) + geom_point(size = 1) +  scale_color_manual(values = c("magenta", "darkgreen", "gold", "dodgerblue", "red")) + labs(title = "PCoA/Bray-Curtis")
dev.off()

#################
### WINDOWS of switch time
#################
# Plot switch right before and right after diet switch
comp = subset_samples(ps, Treatment != "AL_DR20M")
comp = subset_samples(comp, Treatment != "AL_DR12M")
comp = subset_samples(comp, Months == 16)
uf_distance_small  <- phyloseq::distance(comp, "uUniFrac")
wuf_distance_small <- phyloseq::distance(comp, "wUniFrac")
# Save distance matrices for python
write.table(as.data.frame(as.matrix(uf_distance_small)), paste0(path, "analysis/intermediate/16MSW_unifraq_unweighted_dist_", new_day, ".csv"), sep = ";", quote = F)
write.table(as.data.frame(as.matrix(wuf_distance_small)), paste0(path, "analysis/intermediate/16MSWunifraq_weighted_dist_", new_day, ".csv"), sep = ";", quote = F)
# Add meta for analysis
df_d <- melt(as.matrix(uf_distance_small), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
write.table(df_d_new, file = paste0(path, "analysis/intermediate/16MSWunifraq_unweighted-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
df_d <- melt(as.matrix(wuf_distance_small), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
write.table(df_d_new, file = paste0(path, "analysis/intermediate/16MSWunifraq_weighted-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
# Plot
uf_ord <- ordinate(ps, method="PCoA", distance=uf_distance_small)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/16MSWunifraq_unweighted_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, uf_ord, color="Treatment") + theme(aspect.ratio=1, axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=18)) + geom_point(size = 3) + scale_color_manual(values = c("darkgreen", "dodgerblue", "red")) + labs(title = "PCoA/Unweighted UniFraq - 16 Months")
dev.off()

wuf_ord <- ordinate(ps, method="PCoA", distance = wuf_distance_small)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/16MSWunifraq_weighted_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, wuf_ord, color="Treatment") + theme(aspect.ratio=1, axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=18)) + geom_point(size = 3) +  scale_color_manual(values = c("darkgreen", "dodgerblue", "red")) + labs(title = "PCoA/Weighted UniFraq - 16 Months")
dev.off()

#################
### WINDOWS of switch time
#################
# Plot switch right before and right after diet switch
comp = subset_samples(ps, Treatment != "AL_DR16M")
comp = subset_samples(comp, Treatment != "AL_DR12M")
comp = subset_samples(comp, Months == 20)
uf_distance_small  <- phyloseq::distance(comp, "uUniFrac")
wuf_distance_small <- phyloseq::distance(comp, "wUniFrac")
# Save distance matrices for python
write.table(as.data.frame(as.matrix(uf_distance_small)), paste0(path, "analysis/intermediate/20MSW_unifraq_unweighted_dist_", new_day, ".csv"), sep = ";", quote = F)
write.table(as.data.frame(as.matrix(wuf_distance_small)), paste0(path, "analysis/intermediate/20MSWunifraq_weighted_dist_", new_day, ".csv"), sep = ";", quote = F)
# Add meta for analysis
df_d <- melt(as.matrix(uf_distance_small), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
write.table(df_d_new, file = paste0(path, "analysis/intermediate/20MSWunifraq_unweighted-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
df_d <- melt(as.matrix(wuf_distance_small), varnames = c("row", "col"))
df_d_new <- data.frame(df_d)
write.table(df_d_new, file = paste0(path, "analysis/intermediate/20MSWunifraq_weighted-melt_", new_day, ".csv"), col.names = TRUE, sep = ";", quote = F)
# Plot
uf_ord <- ordinate(ps, method="PCoA", distance=uf_distance_small)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/20MSWunifraq_unweighted_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, uf_ord, color="Treatment") + theme(aspect.ratio=1, axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=18)) + geom_point(size = 3) + scale_color_manual(values = c("gold", "dodgerblue", "red")) + labs(title = "PCoA/Unweighted UniFraq - 20 Months")
dev.off()

wuf_ord <- ordinate(ps, method="PCoA", distance = wuf_distance_small)
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/20MSWunifraq_weighted_PCoA_TreatColor_", new_day, ".pdf"))
plot_ordination(ps, wuf_ord, color="Treatment") + theme(aspect.ratio=1, axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=18)) + geom_point(size = 3) +  scale_color_manual(values = c("gold", "dodgerblue", "red")) + labs(title = "PCoA/Weighted UniFraq - 20 Months")
dev.off()

#################
## Checking only AL_DR16M distribution
#################
comp = subset_samples(ps, Treatment == "AL_DR16M")
bray_comp = phyloseq::distance(comp, method="bray", weighted=F)
bray_ord = ordinate(comp, method="PCoA", distance=bray_comp)

# on plot_ordination, to link between collections of same sample
#plot_ordination(comp, bray_ord, "ID", color="Months") + theme(aspect.ratio=1) + geom_path()
library(RColorBrewer)
cols <- rev(brewer.pal(17, 'Greens'))
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/beta_PCoA_ALDR16M_", new_day, ".pdf"))
plot_ordination(comp, bray_ord, color="Months", title = "AL_DR16M") + theme(aspect.ratio=1) +
scale_colour_gradientn(colours = cols)
dev.off()

#################
### Stacked bar plots of taxa abundance
#################
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
# Adding geom_bar(stat="identity") To remove black lines around the boxes of the plots
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/barplot_topOTU_", new_day, ".pdf"))
plot_bar(ps.top20, x="Months", fill="Family") + facet_wrap(~Treatment) + geom_bar(stat="identity")
dev.off()

#################
### Heatmap of just firmicutes
#################
comp = subset_taxa(ps, (Phylum =="Firmicutes") | is.na(Phylum))
# Find which rows (samples) have zero for everything
which(rowSums(otu_table(comp), na.rm = TRUE) == 0)
# Remove them
comp = subset_samples(comp, sample_names(comp) != "Food-7" & sample_names(comp) != "SGRO-0673-5" & sample_names(comp) != "SGRO-0421-4" & sample_names(comp) != "SGRO-0672-5" & sample_names(comp) != "SGRO-0369-12")
#gpt <- prune_taxa(names(sort(taxa_sums(ps), decreasing=TRUE))[1:100], ps)

pdf(paste0(path,"analysis/plots/alpha_beta_diversity/Heatmap_Family_", new_day, ".pdf"))
(p <- plot_heatmap(comp, "NMDS", "bray", "Treatment", "Family", low="#66CCFF", high="#000033", na.value="white"))
dev.off()

#################
### Rarecurves
#################
pdf(paste0(path,"analysis/plots/alpha_beta_diversity/rarecurve_", new_day, ".pdf"))
rarecurve(t(count_tab), step=20, col=sam$Diet_color, lwd=2, ylab="ASVs", label=F) #, xlim = c(0, ))
abline(v=(min(rowSums(t(count_tab)))))
dev.off()
