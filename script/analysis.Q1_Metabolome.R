##########################################################################
##########################################################################
# --------------------------------------------------------------------- #
# ////// Question 1 Script. Native vs nonnative grass metabolomes ///// #
# --------------------------------------------------------------------- #
##########################################################################
##########################################################################
# Colin Richard Morrison
# The University of Texas at Austin, Department of Integrative Biology 
# Brackenridge Field Laboratory, 2907 Lake Austin Blvd., Austin TX 78703
# *crmorrison@utexas.edu# 
# additional analyses for questions 2 and 3 paper by Dr. Liz Bowman and Colin Morrison accompany this script on GitHub 
# updated January 2026

### Goal: Test if native and nonnative grass have distinct metabolomic profiles

### laoad metabolome data 
metab<-load("data/G045_v2_jaccard_total_classified_90_pathway_prob_20251017.RData")

### load host plant metadata
meta<-read.csv("data/G045.Field.data.csv")


################################################################################
### --- 1. PCA  to visualize native v non-native metabolomes ordinations --- ###
################################################################################

### Align metadata to samples x metabolites matrix
# Convert sample IDs to characters (prevents factor/numeric mismatch)
rownames(sampsByCompounds_bin) <- as.character(rownames(sampsByCompounds_bin))
meta$SampleID <- as.character(meta$SampleID)
# Reorder meta to match the rows
meta_ordered <- meta[match(rownames(sampsByCompounds_bin), meta$SampleID), ]
# Create the group vector 
group_var <- meta_ordered$Native.status # change between PC1 and PC2


### ---  PCA on presence/absence metabolite data ###
# Identify columns of metaboites with no variance (metab feature in all samples)
zero_var_cols <- apply(sampsByCompounds_bin, 2, function(x) var(as.numeric(x)) == 0)
sum(zero_var_cols)   # 136/1723 metabolites = present all samples --> omitted to permit centering for PCA
samps_clean <- sampsByCompounds_bin[, !zero_var_cols]

# Hellinger transform matrix
samps_hell <- decostand(samps_clean, method = "hellinger")

# Run PCA
pca_res <- prcomp(samps_hell, scale. = FALSE)
scores <- as.data.frame(pca_res$x)
scores$Group <- group_var
scores$SampleID <- rownames(sampsByCompounds_bin)

# Extract proportion of variance explained
var_explained <- summary(pca_res)$importance[2, ] * 100 


### plot it
# FIGURE 1A: NATIVE vs NON-native metabolomes
pca1<-ggplot(scores, aes(PC1, PC2, color = Group, fill = Group)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.7, alpha = 0.95, color = "black") +
  stat_ellipse(aes(group = Group), level = 0.95, size = 0.8, alpha = 0.5,linetype="longdash") +
  xlab(paste0("PC1 (", format(round(var_explained[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC2 (", format(round(var_explained[2], 1), nsmall = 1), "%)")) +
  scale_color_manual(values = c("Native"= "#00688B",
                                "Nonnative" = "#D3D3D3")) +
  scale_fill_manual(values = c("Native" = "#00688B",
                               "Nonnative" = "#D3D3D3" )) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 20), vjust = 0.5, size = 16),
    axis.title.y = element_text(margin = margin(t = 20), vjust = 5.0, size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    plot.margin  = margin(t = 5, r = 0, b = 10, l = 25, unit = "pt")
  )
#
pca1


### SUPPLEMENTAL FIGURE: NATIVE vs Established vs NON-native 
pca2<-ggplot(scores, aes(PC1, PC2, color = Group, fill = Group)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.7, alpha = 0.95, color = "black") +
  stat_ellipse(aes(group = Group), level = 0.95, size = 0.8, alpha = 0.5,linetype="longdash") +
  xlab(paste0("PC1 (", format(round(var_explained[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC2 (", format(round(var_explained[2], 1), nsmall = 1), "%)")) +
  scale_color_manual(values = c("Native"  = "#00688B",
                                "Established" = "#A52A2A",
                                "Widespread.invasive" = "#D3D3D3")) +
  scale_fill_manual(values = c("Native" = "#00688B",
                               "Established" = "#A52A2A",
                               "Widespread.invasive" = "#D3D3D3" )) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 20), vjust = 0.5, size = 16),
    axis.title.y = element_text(margin = margin(t = 20), vjust = 5.0, size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    plot.margin  = margin(t = 5, r = 0, b = 10, l = 25, unit = "pt")
  )
pca2

ggsave("figure/SupplementaryFigureS2.jpeg",
       plot = pca2,
       width = 12, height = 6,
       device = 'jpeg',
       dpi = 600)


################################################################################################
### --- 2. Assess grass metabolome variation as function of native status with PERMANOVA --- ###
################################################################################################

rownames(jaccard) <- as.character(rownames(jaccard))
meta$SampleID <- as.character(meta$SampleID)
# Reorder meta to match the rows 
meta_ordered <- meta[match(rownames(jaccard), meta$SampleID), ]
# Create the group vector
group_var <- meta_ordered$Native.status.fine # change for different nativeness category


### Run PERMANOVA on NATIVE v NON-native
permanova_res <- adonis2(jaccard ~ Native.status * Genus,
                         by = "term",
                         data = meta_ordered, 
                         permutations = 9999)
permanova_res


### Run pairwise PERMANOVAs, with hom-bonferroni correct, comparing NATIVE v Non-native v Established
pairwise_permanova <- function(jaccard, group_var, permutations = 9999, method = "holm") {
  groups <- unique(group_var)
  results <- data.frame()
  for(i in 1:(length(groups)-1)) {
    for(j in (i+1):length(groups)) {
      pair <- groups[c(i,j)]
      idx <- which(group_var %in% pair)
      sub_dist <- as.dist(as.matrix(jaccard)[idx, idx])
      sub_group <- group_var[idx]
      temp_meta <- data.frame(sub_group = sub_group)
      res <- adonis2(sub_dist ~ sub_group, data = temp_meta, permutations = permutations)
      # Add group row
      results <- rbind(results,
                       data.frame(
                         Pair = paste(pair, collapse = " vs "),
                         Term = rownames(res)[1],
                         Df = res$Df[1],
                         SumOfSqs = res$SumOfSqs[1],
                         R2 = res$R2[1],
                         F = res$F[1],
                         p = res$`Pr(>F)`[1]
                       ))
      # Add residual row
      results <- rbind(results,
                       data.frame(
                         Pair = paste(pair, collapse = " vs "),
                         Term = rownames(res)[2],
                         Df = res$Df[2],
                         SumOfSqs = res$SumOfSqs[2],
                         R2 = res$R2[2],
                         F = NA,
                         p = NA
                       ))
    }
  }
  # Apply holm-bonferroni p-value adjustment on group rows
  group_idx <- !is.na(results$p)
  results$p_adj <- NA
  results$p_adj[group_idx] <- p.adjust(results$p[group_idx], method = method)
  return(results)
}

# Run pairwise PERMANOVA on your Jaccard matrix
pairwise_res2 <- pairwise_permanova(jaccard, group_var, permutations = 9999)
pairwise_res2



##################################################################################
### --- 3.Calculate metabolomic richness of each grass sample and species --- ###
##################################################################################

### metabolomic richness of each sample
sample_metab_rich <- sampsByCompounds_bin %>% 
  as.data.frame() %>%                   
  rownames_to_column("SampleID") %>%    
  mutate(metabolite_richness = rowSums(across(-SampleID))) %>% 
  select(SampleID, metabolite_richness)
sample_metab_rich


### average metabolomic richness of each grass species +/- standard error
meta$SampleID <- as.character(meta$SampleID)
richness_summary$SampleID <- as.character(richness_summary$SampleID)

meta_richness <- meta %>%
  left_join(sample_metab_rich, by = "SampleID")

species_metab_rich <- meta_richness %>%
  group_by(Species_final_nov24) %>%
  summarise(
    n_samples = n(),
    mean_richness = mean(metabolite_richness.x, na.rm = TRUE),
    sd_richness = sd(metabolite_richness.x, na.rm = TRUE),
    se_richness = sd_richness / sqrt(n_samples)
  )
species_metab_rich

### - END - ##############################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################