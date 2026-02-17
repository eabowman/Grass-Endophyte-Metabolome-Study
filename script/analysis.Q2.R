# Assessment of data to answer Q2. How does the endophyte community correlate
# to metabolomics?
# Written by Liz Bowman, PhD eabowman@utexas.edu
# Dec. 19, 2024

# Host data ----
# Phylogenetic tree
host.tree <- read.tree("data/G045_Concatenated_ITStrnL_withbranchlengths.nwk")

# Calculate host phylogenetic distance
dist.host <- cophenetic.phylo(host.tree)
dim(dist.host)

# Field data
data.host.sp <- read.csv('data/G045.Field.data.csv', header = T)

## Compute PCNM eigenvectors for location data----
# Geographic distance
dist.geo <- vegdist(data.host.sp[c('lat','long')], method = 'euclidean')
geo.pcnm <- pcnm(dist.geo)

# add to us.clim data frame
data.host.sp$geo.pcnm.1 <- geo.pcnm$vectors[,1]

# Culture data ----
data.cult <- read.csv('data/G045.CultureBased.data.csv')

# OTU data, culture-free data ----
data.otu <- read.csv('data/G045.CultureFree.data_NegRemoved_rarefied.csv')

# make SampleID column in data.otu an integer
data.otu$SampleID <- as.character(data.otu$SampleID)

# isolate endophyte community data 
data.comm <- data.otu[-1]

# Remove columns with 0 or singletons
data.comm[colSums(data.comm) > 1] -> data.comm

# Calculate fungal fisher's alpha from culture free data
data.cult %>%
  mutate(read.fisher.alpha = fisher.alpha(data.comm),
         log.read.fisher.alpha = log(fisher.alpha(data.comm)),
         read.abundance <- rowSums(data.comm),
         read.richness <- specnumber(data.comm)) -> data.all

# Metabolomic distance table ----
## Overall metabolomics ----
data.metab <- read.csv('data/G045.Metabolomic.data_reorganized_JaccardV2.csv',
                       row.names = 'SampleID')

# Clean up distance table
colnames(data.metab) <- rownames(data.metab)
data.metab %>%
  mutate(SampleID = as.character(rownames(.))) -> data.metab

## Terpinoids -----
load(file = 'data/G045_v2_jaccard_terpenoids_20251017.RData')
data.terp <- jaccard
colnames(data.terp) <- rownames(data.terp)
data.terp %>%
  mutate(SampleID = as.character(rownames(.))) -> data.terp

## Shikimates -----
load(file = 'data/G045_v2_jaccard_shikimates_20251017.RData')
data.shik <- jaccard
colnames(data.shik) <- rownames(data.shik)
data.shik %>%
  mutate(SampleID = as.character(rownames(.))) -> data.shik

## Fatty Acids -----
load(file = 'data/G045_v2_jaccard_fatty_acids_20251023.RData')
data.fat <- jaccard
colnames(data.fat) <- rownames(data.fat)
data.fat %>%
  mutate(SampleID = as.character(rownames(.))) -> data.fat

# 1. How does the endophytes community correlate to the metabolomics? ----
# To assess the correlation of the endophyte community to the metabolome a 
# distance based redundancy analysis (db-RDA). 
# https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/

# Calculate culture-free distance matrices with occurrence- (Jaccard) and
# abundance-based (Morisita Horn) measures

# Hellinger transformation
comm.hell <- decostand(data.comm, method = 'total')

# Jaccard dissimilarity
dist.jacc <- vegdist(comm.hell, method = 'jaccard', upper = F) 

# Morisita-Horn dissimilarity
dist.horn <- vegdist(comm.hell, method = 'horn', upper = F) 

# convert data.metab to a distance matrix
data.metab %>% select(-SampleID) -> data.metab.dist
as.dist(data.metab.dist, upper = F) -> dist.metab
rm(data.metab.dist)

# Turn dist.host into a distance table with class 'dist'
dist(dist.host) -> dist.host

# covert geo.pcnm.1 to a distance matrix
vegdist(data.host.sp$geo.pcnm.1,
        method = 'euclidean') -> dist.geo

# covert Native.status to a distance matrix
data.host.sp %>%
  mutate(Native.status.bin = 
           case_when(Native.status == 'Native' ~ 0,
                     Native.status == 'Nonnative' ~ 1)) -> data.host.sp
vegdist(data.host.sp$Native.status.bin,
        method = 'euclidean',
        binary = T) -> dist.native

## MRM ----
mrm.df <- ecodist::MRM(formula = dist.metab ~  dist.host + dist.horn + dist.native,
             nperm = 1000, method = "linear", mrank = T)

as.data.frame(mrm.df$coef) -> mrm.coef
write.csv(mrm.coef,
          'results/MultipleRegressionDistanceMatrix_Coefficient.csv', 
          row.names = T)

as.data.frame(mrm.df$r.squared) -> mrm.r2
write.csv(mrm.r2,
          'results/MultipleRegressionDistanceMatrix_Rsquared.csv', 
          row.names = T)

as.data.frame(mrm.df$F.test) -> mrm.Ftest
write.csv(mrm.Ftest,
          'results/MultipleRegressionDistanceMatrix_Ftest.csv', 
          row.names = T)


## Partial Mantel test ----
### Overall ----
### how much does host genetic distance contribute to metabolomic diversity
### when variation due to the endophyte community is  removed?

all.host.horn <- mpmcorrelogram(dist.metab, dist.host, dist.horn,
                               method = 'spearman')

# Plot in ggplot
all.host.horn.df <- data.frame(class = seq_along(all.host.horn$clases),
                             r = all.host.horn$rM,
                             p = all.host.horn$pval.Bonferroni,
                             dist = all.host.horn$breaks[-length(all.host.horn$breaks)])
mutate(all.host.horn.df, sig = p < 0.05) -> all.host.horn.df

Fig.2c <- ggplot(all.host.horn.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

### how much does the endophyte community contribute to metabolomic diversity
### when variation due to host genetic diversity is removed?

all.horn.host <- mpmcorrelogram(
  dist.metab, dist.horn, dist.host,
  nclass = 6, method = "spearman"
)

# Plot in ggplot
all.horn.host.df <- data.frame(class = seq_along(all.horn.host$clases),
                             r = all.horn.host$rM,
                             p = all.horn.host$pval.Bonferroni,
                             dist = all.horn.host$breaks[-length(all.horn.host$breaks)])
mutate(all.horn.host.df, sig = p < 0.05) -> all.horn.host.df

Fig.2d <- ggplot(all.horn.host.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

### Subsets of metabolites ----

### data.terp 
# convert data.terp to a distance matrix
data.terp %>% select(-SampleID) -> data.terp.dist
as.dist(data.terp.dist, upper = F) -> dist.terp
rm(data.terp.dist)

### how much does the endophyte community contribute to metabolomic diversity
### when variation due to host genetic diversity is removed?

terp.horn.host <- mpmcorrelogram(
  dist.terp, dist.horn, dist.host,
  nclass = 6, 
  method = "spearman"
)

# Plot in ggplot
terp.horn.host.df <- data.frame(class = seq_along(terp.horn.host$clases),
                             r = terp.horn.host$rM,
                             p = terp.horn.host$pval.Bonferroni,
                             dist = terp.horn.host$breaks[-length(terp.horn.host$breaks)])
mutate(terp.horn.host.df, sig = p < 0.05) -> terp.horn.host.df

Supplfig.6a <- ggplot(terp.horn.host.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

### how much does host genetic distance contribute to metabolomic diversity
### when variation due to the endophyte community is  removed?

terp.host.horn <- mpmcorrelogram(dist.terp, dist.host, dist.horn,
                               method = 'spearman')

# Plot in ggplot
terp.host.horn.df <- data.frame(class = seq_along(terp.host.horn$clases),
                             r = terp.host.horn$rM,
                             p = terp.host.horn$pval.Bonferroni,
                             dist = terp.host.horn$breaks[-length(terp.host.horn$breaks)])
mutate(terp.host.horn.df, sig = p < 0.05) -> terp.host.horn.df

Supplfig.6b <- ggplot(terp.host.horn.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

### data.shik
# convert data.shik to a distance matrix
data.shik %>% select(-SampleID) -> data.shik.dist
as.dist(data.shik.dist, upper = F) -> dist.shik
rm(data.shik.dist)

### how much does the endophyte community contribute to metabolomic diversity
### when variation due to host genetic diversity is removed?

Shik.horn.host <- mpmcorrelogram(
  dist.shik, dist.horn, dist.host,
  nclass = 6, 
  method = "spearman"
)

# Plot in ggplot
Shik.horn.host.df <- data.frame(class = seq_along(Shik.horn.host$clases),
                             r = Shik.horn.host$rM,
                             p = Shik.horn.host$pval.Bonferroni,
                             dist = Shik.horn.host$breaks[-length(Shik.horn.host$breaks)])
mutate(Shik.horn.host.df, sig = p < 0.05) -> Shik.horn.host.df

Supplfig.6c <- ggplot(Shik.horn.host.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

### how much does host genetic distance contribute to metabolomic diversity
### when variation due to the endophyte community is  removed?

shik.host.horn <- mpmcorrelogram(dist.shik, dist.host, dist.horn,
                               method = 'spearman')

# Plot in ggplot
shik.host.horn.df <- data.frame(class = seq_along(shik.host.horn$clases),
                             r = shik.host.horn$rM,
                             p = shik.host.horn$pval.Bonferroni,
                             dist = shik.host.horn$breaks[-length(shik.host.horn$breaks)])
mutate(shik.host.horn.df, sig = p < 0.05) -> shik.host.horn.df

Supplfig.6d <- ggplot(shik.host.horn.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

### data.fat
# convert data.fat to a distance matrix
data.fat %>% select(-SampleID) -> data.fat.dist
as.dist(data.fat.dist, upper = F) -> dist.fat
rm(data.fat.dist)

### how much does the endophyte community contribute to metabolomic diversity
### when variation due to host genetic diversity is removed?

fat.horn.host <- mpmcorrelogram(
  dist.fat, dist.horn, dist.host,
  nclass = 6, 
  method = "spearman"
)

# Plot in ggplot
fat.horn.host.df <- data.frame(class = seq_along(fat.horn.host$clases),
                             r = fat.horn.host$rM,
                             p = fat.horn.host$pval.Bonferroni,
                             dist = fat.horn.host$breaks[-length(fat.horn.host$breaks)])
mutate(fat.horn.host.df, sig = p < 0.05) -> fat.horn.host.df

Supplfig.6e <- ggplot(fat.horn.host.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')


### how much does host genetic distance contribute to metabolomic diversity
### when variation due to the endophyte community is  removed?

fat.host.horn <- mpmcorrelogram(dist.fat, dist.host, dist.horn,
                               method = 'spearman')

# Plot in ggplot
fat.host.horn.df <- data.frame(class = seq_along(fat.host.horn$clases),
                             r = fat.host.horn$rM,
                             p = fat.host.horn$pval.Bonferroni,
                             dist = fat.host.horn$breaks[-length(fat.host.horn$breaks)])
mutate(fat.host.horn.df, sig = p < 0.05) -> fat.host.horn.df

Supplfig.6f <- ggplot(fat.host.horn.df,
                              aes(x = dist,
                                  y = r)) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_segment(aes(xend = dist,
                   yend = 0)) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Distance class",
       y = expression(r[M])) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none')

# Generate Supplementary Fig. S6
Supplfig.6a +
  Supplfig.6b +
  Supplfig.6c +
  Supplfig.6d +
  Supplfig.6e +
  Supplfig.6f +
    plot_layout(
    nrow = 3,
    ncol = 2
  ) + 
  plot_annotation(tag_levels = "A") -> SupplFig.S6

ggsave("figure/SupplementaryFigureS6.jpeg",
       plot = SupplFig.S6,
       width = 10, height = 10,
       device = 'jpeg',
       dpi = 600)

## VarPart: Table 1 ----

# Convert distances to PCoA axes
metab.pcoa <- cmdscale(dist.metab, k = 5, eig = TRUE)$points  # response
endo.pcoa  <- cmdscale(dist.horn,  k = 5, eig = TRUE)$points  # predictor 1
host.pcoa  <- cmdscale(dist.host,  k = 5, eig = TRUE)$points  # predictor 2

var.metab <- varpart(metab.pcoa, endo.pcoa, host.pcoa)
plot(var.metab)

# Unique contribution of endophytes controlling for host
rda.endo <- rda(metab.pcoa, endo.pcoa, host.pcoa)
anova(rda.endo, by = "axis")

# Unique contribution of host controlling for endophytes
rda.host <- rda(metab.pcoa, host.pcoa, endo.pcoa)
anova(rda.host, by = "axis")

## Subsets: Need to run partial mantel section (above) for subsets to run this section
# Convert distances to PCoA axes
terp.pcoa <- cmdscale(dist.terp, k = 5, eig = TRUE)$points  # response
var.terp <- varpart(terp.pcoa, endo.pcoa, host.pcoa)
plot(var.terp) #34% FE, 25% host

shik.pcoa <- cmdscale(dist.shik, k = 5, eig = TRUE)$points  # response
var.shik <- varpart(shik.pcoa, endo.pcoa, host.pcoa)
plot(var.shik) #24% FE, 28% host

fat.pcoa <- cmdscale(dist.fat, k = 5, eig = TRUE)$points  # response
var.fat <- varpart(fat.pcoa, endo.pcoa, host.pcoa)
plot(var.fat) #20% FE, 26% host

## Linear regression plots ----
### organize data
#### Morisita-Horn
as.data.frame(as.matrix(dist.horn)) -> df.horn.dist
colnames(df.horn.dist) <- data.otu$SampleID
df.horn.dist$sample1 <- data.otu$SampleID

df.horn.dist %>%
  gather(key = sample2, value = dist.horn, -sample1) %>%
  filter(sample1 > sample2) %>%
  mutate(sample1 = as.character(sample1)) -> df.horn.dist

#### Jaccard
as.data.frame(as.matrix(dist.jacc)) -> df.jacc.dist
colnames(df.jacc.dist) <- data.otu$SampleID
df.jacc.dist$sample1 <- data.otu$SampleID

df.jacc.dist %>%
  gather(key = sample2, value = dist.jacc, -sample1) %>%
  filter(sample1 > sample2) %>%
  mutate(sample1 = as.character(sample1)) -> df.jacc.dist

#### host genetic distance
as.data.frame(as.matrix(dist.host)) -> dist.host
dist.host$sample1 <- colnames(dist.host)

dist.host %>%
  gather(key = sample2, value = dist.host, -sample1) %>%
  mutate(sample2 = stringr::str_remove(sample2, "^X"),
         sample1 = as.character(sample1)) %>%
  filter(sample1 > sample2) -> dist.host.long

#### metabolomic distance
as.data.frame(as.matrix(dist.metab)) -> dist.metab
dist.metab$sample1 <- colnames(dist.metab)

dist.metab %>%
  gather(key = sample2, value = dist.metab, -sample1) %>%
  mutate(sample2 = stringr::str_remove(sample2, "^X"),
         sample1 = as.character(sample1)) %>%
  filter(sample1 > sample2) -> dist.metab.long

#### Join all data fdist.host.gd.long#### Join all data frames ----
full_join(df.horn.dist, df.jacc.dist) -> dist.all

full_join(dist.all, dist.host.long) -> dist.all
full_join(dist.all, dist.metab.long) -> dist.all

#### Add groups ----
# Native status
for(i in unique(data.all$SampleID)){
  Nat.stat.i <- data.all[data.all$SampleID == i, 'Native.status']
  dist.all[dist.all$sample1 == i, "Native.Status.sample1"] <- Nat.stat.i
  dist.all[dist.all$sample2 == i, "Native.Status.sample2"] <- Nat.stat.i
}
dist.all$Native.Status.group <- NA
for(i in 1:nrow(dist.all)){
  if(dist.all[i, 'Native.Status.sample1'] == dist.all[i, 'Native.Status.sample2']){
    dist.all[i, 'Native.Status.group'] <- 'Within'
  } else if(dist.all[i, 'Native.Status.sample1'] != dist.all[i,
                                                             'Native.Status.sample2']){
    dist.all[i, 'Native.Status.group'] <- 'Between'
  }}

# Genus
for(i in unique(data.all$SampleID)){
  Gen.stat.i <- data.all[data.all$SampleID == i, 'Genus']
  dist.all[dist.all$sample1 == i, "Genus.sample1"] <- Gen.stat.i
  dist.all[dist.all$sample2 == i, "Genus.sample2"] <- Gen.stat.i
}
dist.all$Genus.group <- NA
for(i in 1:nrow(dist.all)){
  if(dist.all[i, 'Genus.sample1'] == dist.all[i, 'Genus.sample2']){
    dist.all[i, 'Genus.group'] <- 'Within'
  } else if(dist.all[i, 'Genus.sample1'] != dist.all[i,
                                                             'Genus.sample2']){
    dist.all[i, 'Genus.group'] <- 'Between'
}}

# Check normality
hist(dist.all$dist.horn)
hist(dist.all$dist.jacc)
hist(dist.all$dist.metab)
shapiro.test(dist.all$dist.metab)
hist(dist.all$dist.host)
shapiro.test(log1p(dist.all$dist.host))

#### Metabolome as a function of the Endophyte community -----
as.factor(dist.all$Native.Status.group) -> dist.all$Native.Status.group

Fig.2a <- dist.all %>%
  filter(Native.Status.group == 'Within') %>%
  ggplot(aes(x = dist.horn,
           y = dist.metab,
           fill = dist.host)) +
  geom_point(shape = 21,
             color = 'black',
             size = 1.5) +
  # facet_grid(Native.Status.sample1 ~ Native.Status.sample2) +
  # scale_fill_manual(values = c('#086c8d', '#d2d2d2')) +
  geom_smooth(method = 'lm',
              alpha = 0.6,
              color = '#f39b9b',
              linewidth = 0.5) +
  scale_fill_gradient(low = '#abdbe3', high = '#063970') +
  xlab('Endophyte community dissimilarity') +
  ylab('Metabolic dissimilarity') +
  ylim(0.1, 0.8) +
  labs(fill = 'Host\nphylogenetic\ndistance') +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'right',
        legend.title = element_text(size = 10))

#### Metabolome as a function of host plant -----
Fig.2b <- dist.all %>%
  ggplot(aes(x = dist.host,
           y = dist.metab))+
  geom_point(shape = 21,
             color = 'black',
             fill = 'darkgrey',
             size = 1.5) +
  # scale_fill_manual(values = c('#086c8d', '#d2d2d2')) +
  geom_smooth(method = 'lm',
              alpha = 0.6,
              color = '#f39b9b',
              linewidth = 0.5) +
  facet_grid(Native.Status.sample1 ~ Native.Status.sample2) +
  xlab('Host phylogenetic dissimilarity') +
  ylab('Metabolomic distance') +
  ylim(0.1, 0.8) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none',
        strip.text = element_text(size = 12, color = 'black'),
        strip.background = element_blank())

# Make Figure 2
Fig.2a +
  Fig.2b +
  Fig.2c +
  Fig.2d +
    plot_layout(
    ncol = 2,
    nrow = 2,
    byrow = T
  ) + 
  plot_annotation(tag_levels = "A") -> Fig.2

ggsave("figure/Figure2.jpeg",
       plot = Fig.2,
       width = 12, height = 6,
       device = 'jpeg',
       dpi = 600)

