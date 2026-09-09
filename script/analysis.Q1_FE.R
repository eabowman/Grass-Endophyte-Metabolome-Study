# Analyses to answer Q1: How does the endophyte community vary by host plant?
# Written by Liz Bowman, PhD eabowman@utexas.edu
# Sep. 20, 2024

# Metadata ----
data.host.sp <- read.csv('data/G045.Field.data.csv', header = T)

# OTU data, culture-free data ----

## OTU data with PCR and extraction controls removed and rarefied ----
data.free <- read.csv('data/G045.CultureFree.data_NegRemoved_rarefied.csv')

# make SampleID column in data.otu an integer
data.free$SampleID <- as.integer(data.free$SampleID)

data.host.sp %>%
  select(SampleID, Genus, Species_final_nov24, Native.status, Native.status.fine,
         lat,long) %>%
  full_join(., data.free, by = 'SampleID') -> data.free

# 1. How does the endophyte community vary by host plant? ----
# Here, we are considering both native status and host genetic distance

## Culture-free ----
# isolate endophyte community data 
data.free.comm <- select(data.free, starts_with('Otu'))

# Remove columns with 0 or singletons
data.free.comm[colSums(data.free.comm) > 1] -> data.free.comm

### distance based Moran's eigenvector (dbMEM) -----
# Load varpart2.MEM.R function file
# varpart2.MEM.R was written by Legendre et al. 2012 'Variation partitioning
# involving orthogonal spatial eigenfunction submodels'
source('script/varpart2.MEM.R')

# isolate location data
gps.xy <- data.free[c('lat','long')]

### PCNM
# Create distance matrix of spatial data
gps.dist <- dist(gps.xy, method = 'euclidean')

# transform spatial data 
gps.pcnm <- pcnm(gps.dist)

# access eigenvectors with scores() function

# map spatial data
op <- par(mfrow = c(1,2))
ordisurf(data.free, scores(gps.pcnm, choices=1), bubble = 4, main = "PCNM 1")
ordisurf(data.free, scores(gps.pcnm, choices=2), bubble = 4, main = "PCNM 2")
par(op)

data.free$PCNM1 <- gps.pcnm$vectors[,1]
data.free$PCNM2 <- gps.pcnm$vectors[,2]

# Culture based data ----
# read in data frame created below; all samples (nrows = 36)
data.cult <- read.csv('data/G045.CultureBased.data.csv')

# Isolate only rows with culture based sequence data (nrows = 18)
data.cult.short <- filter(data.cult, !is.na(cOtu0))

# isolate culture community
data.cult.short %>%
  select(starts_with('cOtu')) -> data.cult.comm

### distance based Moran's eigenvector (dbMEM) -----

# isolate location data
gps.xy <- data.cult.short[c('lat','long')]

### PCNM
# Create distance matrix of spatial data
gps.dist <- dist(gps.xy, method = 'euclidean')

# transform spatial data 
gps.pcnm <- pcnm(gps.dist)

# access eigenvectors with scores() function

# map spatial data
op <- par(mfrow = c(1,2))
ordisurf(data.free, scores(gps.pcnm, choices=1), bubble = 4, main = "PCNM 1")
ordisurf(data.free, scores(gps.pcnm, choices=2), bubble = 4, main = "PCNM 1")
par(op)

data.cult.short$PCNM1 <- gps.pcnm$vectors[,1]
data.cult.short$PCNM2 <- gps.pcnm$vectors[,2]

# Diversity and abundance ----
# calculate diversity of culture-free endophyte community: shannon
data.free$read.Shannon <- diversity(data.free.comm, index = 'shannon')

# calculate diversity of culture-free endophyte community: Fisher's alpha
data.free$read.Fishers.alpha <- fisher.alpha(data.free.comm)
data.free$log.read.Fishers.alpha <- log(data.free$read.Fishers.alpha)

# Read abundance and OTU richness
CF.abundance <- data.frame(SampleID = data.free$SampleID,
                           read.abundance = rowSums(data.free.comm),
                           read.richness = specnumber(data.free.comm))
full_join(data.free, CF.abundance) -> data.free

# Calculate diversity of culture-based endophyte community: shannon
data.cult.short$cult.shannon <- diversity(data.cult.comm, index = 'shannon')
hist(data.cult.short$cult.shannon)
shapiro.test(data.cult.short$cult.shannon)

# calculate diversity of culture-based endophyte community: Fisher's alpha
data.cult.short$cult.Fishers.alpha <- fisher.alpha(data.cult.comm)

data.cult.short$cult.richness <- specnumber(data.cult.comm)

## Combine culture-based and culture-free data ----
data.cult.short %>%
  mutate(isol.freq = slants.with.growth/total.slants) %>%
  select(SampleID, isol.freq,
         cult.shannon, cult.Fishers.alpha, cult.richness) %>%
  full_join(., data.free, by = 'SampleID') -> data.cf.cb

## Culture-free ----
# OTU richness
hist(data.cf.cb$log.read.Fishers.alpha)
shapiro.test(data.cf.cb$log.read.Fishers.alpha)

cf.aov <- aov(read.richness ~ Native.status + Genus, data = data.cf.cb)
summary(cf.aov)
plot(cf.aov$residuals)

data.cf.cb$Genus <- factor(data.cf.cb$Genus, 
       levels = c('Aristida','Elionurus','Eragrostis','Paspalum','Setaria',
       'Bothriochloa','Cenchrus','Heteropogon','Megathyrsus','Melinis','Urochloa'))

ggplot(data.cf.cb,
       aes(y = read.richness,
           x = Genus,
           fill = Genus,
           color = Native.status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(aes(shape = Native.status),
             size = 2.5) +
  ylab("CF OTU richness") +
  xlab('') +
  scale_fill_brewer(type = 'div', palette = 1) +
  scale_color_manual(values = c('black', 'grey')) +
  scale_shape_manual(values = c(23, 21)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        # axis.text.x = element_text(size = 10, color = 'black',
        #                            angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none') -> cf.div.plot

# Read abundance
hist(data.cf.cb$read.abundance)
shapiro.test(log1p(data.cf.cb$read.abundance))

cf.aov <- aov(read.abundance ~ Native.status + Genus, data = data.cf.cb)
summary(cf.aov)
plot(cf.aov$residuals)

ggplot(data.cf.cb,
       aes(y = read.abundance,
           x = Genus,
           fill = Genus,
           color = Native.status)) +
  geom_point(aes(shape = Native.status),
             size = 2.5) +
  geom_boxplot(alpha = 0.5) +
  ylab("CF read abundance") +
  xlab('') +
  scale_fill_brewer(type = 'div', palette = 1) +
  scale_color_manual(values = c('black', 'grey')) +
  scale_shape_manual(values = c(23, 21)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        # axis.text.x = element_text(size = 10, color = 'black',
        #                            angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12, color = 'black')) -> cf.ab.plot

## Culture-based ----
cb.aov <- aov(cult.richness ~ Native.status + Genus, data = data.cf.cb)
summary(cb.aov)
plot(cb.aov$residuals)

ggplot(data.cf.cb,
       aes(y = cult.richness,
           x = Genus,
           fill = Genus,
           color = Native.status)) +
  geom_point(aes(shape = Native.status),
             size = 2.5) +
  geom_boxplot(alpha = 0.5) +
  ylab("CB OTU richness") +
  xlab('') +
  scale_fill_brewer(type = 'div', palette = 1) +
  scale_color_manual(values = c('black', 'grey')) +
  scale_shape_manual(values = c(23, 21)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none') -> cb.div.plot

# Read abundance
hist(data.cf.cb$isol.freq)
shapiro.test(data.cf.cb$isol.freq)

cb.aov <- aov(isol.freq ~ Native.status + Genus, data = data.cf.cb)
summary(cb.aov)
plot(cb.aov$residuals)

ggplot(data.cf.cb,
       aes(y = isol.freq,
           x = Genus,
           fill = Genus,
           color = Native.status)) +
  geom_point(aes(shape = Native.status),
             size = 2.5) +
  geom_boxplot(alpha = 0.5) +
  ylab("CB isol. frequency") +
  xlab('') +
  scale_fill_brewer(type = 'div', palette = 1) +
  scale_color_manual(values = c('black', 'grey')) +
  scale_shape_manual(values = c(23, 21)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 12, color = 'black'),
        legend.position = 'none') -> cb.ab.plot

# Supplementary Figure S3 ----

(cf.div.plot + cf.ab.plot + cb.div.plot + cb.ab.plot) +
  plot_layout(nrow = 2, ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A") -> SupplFig.S3

ggsave("figure/SupplementaryFigureS3.jpeg",
       plot = SupplFig.S3,
       width = 10, height = 6,
       device = 'jpeg',
       dpi = 600)

# Metabolomic richness as a function of FE abundance and diversity ----
## Culture-free ----
data.free %>%
  select(SampleID, 
         read.Fishers.alpha, log.read.Fishers.alpha, 
         read.abundance, read.richness) %>%
  right_join(data.host.sp) %>%
  mutate(log1p.read.richness = log1p(read.richness)) -> data.host.free

hist(data.host.free$metabolite_richness)
shapiro.test(data.host.free$metabolite_richness)
hist(log1p(data.host.free$read.richness))
shapiro.test(log1p(data.host.free$read.richness))

free.lm <- lm(metabolite_richness ~ log.read.Fishers.alpha,
              data = data.host.free)
summary(free.lm)

ggplot(data.host.free, 
       aes(x = log.read.Fishers.alpha,
           y = metabolite_richness,
           color = Native.status.fine
           )) +
  geom_point() +
  geom_smooth(method = 'glm') +
  theme_classic()

## Culture-based ----
data.host.sp %>%
  filter(SampleID %in% data.cult.short$SampleID) -> data.host.sp.short

data.cult.short %>%
  select(SampleID, 
         isolation.freq, cult.shannon, cult.Fishers.alpha, cult.richness) %>%
    right_join(data.host.sp.short) -> data.host.cult

hist(data.host.cult$metabolite_richness)
shapiro.test(data.host.cult$metabolite_richness)
hist(log1p(data.host.cult$isolation.freq))
shapiro.test(data.host.cult$isolation.freq)


cult.lm <- lm(metabolite_richness ~ cult.richness * Native.status,
              data = data.host.cult)
summary(cult.lm)

ggplot(data.host.cult, 
       aes(x = cult.richness,
           y = metabolite_richness,
           color = Native.status)) +
  geom_point() +
  geom_smooth(method = 'glm') +
  theme_classic()

# PERMANOVA: Jaccard Culture-based ----

# Hellinger transformation
comm.hell <- decostand(data.cult.comm, method = 'total')

# Jaccard dissimilarity
jacc.cult.dist <- vegdist(comm.hell, method = 'jaccard', upper = F,
                          binary = T) 

# Morisita-Horn dissimilarity
horn.cult.dist <- vegdist(comm.hell, method = 'horn', upper = F)

# check for homogeneity of group dispersions
jacc.betadisper <- betadisper(jacc.cult.dist, group = data.cult.short$Genus)
permutest(jacc.betadisper, permutations = 99)
anova(jacc.betadisper)
plot(jacc.betadisper)

jacc.adonis <- adonis2(jacc.cult.dist ~ Native.status + Genus + PCNM1,
                       by = 'term',
                       data = data.cult.short,
                       strata = data.cult.short$site)
jacc.adonis

# write out results
as.data.frame(jacc.adonis) -> jacc.adonis
write.csv(jacc.adonis, 'results/SupplementaryTableS3_CB_Jaccard.csv',
          row.names = T)

#### NMDS plot: Jaccard Culture-based ----
jacc.mds <- metaMDS(jacc.cult.dist, dist = 'bray',
                    try = 1000, trymax = 1000,
                    tidy = T, shrink = T, k = 2)
jacc.stress <- jacc.mds$stress

#### Make data frame for plot ----
data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                          NMDS2 = jacc.mds$points[,2],
                          sample.no = data.cult.short$Sample.no,
                          site = data.cult.short$site,
                          Native.status = data.cult.short$Native.status,
                          Genus = data.cult.short$Genus, 
                          Species = data.cult.short$Species)

jacc.plot.s4b <- ggplot(data = data.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(aes(fill = Genus),
                 shape = 21,
                 color = 'black',
                 size = 2) +  # Using 'Species' for shape too, or adjust as needed
  stat_ellipse(aes(group = Native.status,
                   color = Native.status),  # map color to Native.status
               geom = "polygon",
               linetype = "dashed",
               fill = NA,
               linewidth = 0.8) +
  scale_fill_manual(values = c('#bf812e', '#dfc17e', '#cfeae4', '#81ccc1',
                               '#35968d', '#00665e', '#003b2f')) +
  scale_color_manual(values = c('darkgrey', 'black')) +    
  # geom_text(aes(label = sample.no),  # Change to the desired label column
  #           vjust = -1,              # Adjust vertical position
  #           size = 3) +              # Adjust text size
  coord_equal() +
  labs(color = 'Native status') +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = 'Black'),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = 'none')

# PERMANOVA: Morisita Horn Culture-based ----
# check for homogeneity of group dispersions
horn.betadisper <- betadisper(horn.cult.dist, group = data.cult.short$Genus)
permutest(horn.betadisper, permutations = 99)
anova(horn.betadisper)
plot(horn.betadisper)

horn.adonis <- adonis2(horn.cult.dist ~ Native.status * Genus + PCNM1,
                       by = 'term',
                       data = data.cult.short,
                       strata = data.cult.short$site)
horn.adonis

# write out results
as.data.frame(horn.adonis) -> horn.adonis
write.csv(horn.adonis, 'results/SupplementaryTableS3_CB_MorisitaHorn.csv', row.names = T)

#### NMDS plot: Morisa Horn Culture-based ----
horn.mds <- metaMDS(horn.cult.dist, dist = 'bray',
                    try = 1000, trymax = 1000,
                    tidy = T, shrink = T, k = 2)
horn.stress <- horn.mds$stress

#### Make data frame for plot ----
data.scores <- data.frame(NMDS1 = horn.mds$points[,1],
                          NMDS2 = horn.mds$points[,2],
                          sample.no = data.cult.short$Sample.no,
                          site = data.cult.short$site,
                          Native.status = data.cult.short$Native.status,
                          Genus = data.cult.short$Genus, 
                          Species = data.cult.short$Species)

data.scores$Species <- as.factor(data.scores$Species)

horn.plot.1c <- ggplot(data = data.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(aes(fill = Genus),
                 shape = 21,
                 color = 'black',
                 size = 2) +  # Using 'Species' for shape too, or adjust as needed
  stat_ellipse(aes(group = Native.status,
                   color = Native.status),  # map color to Native.status
               geom = "polygon",
               linetype = "dashed",
               fill = NA,
               linewidth = 0.8) +
  scale_fill_manual(values = c('#bf812e', '#dfc17e', '#cfeae4', '#81ccc1',
                               '#35968d', '#00665e', '#003b2f')) +
  scale_color_manual(values = c('darkgrey', 'black')) +    
  # geom_text(aes(label = sample.no),  # Change to the desired label column
  #           vjust = -1,              # Adjust vertical position
  #           size = 3) +              # Adjust text size
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = 'Black'),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = 'none')

# PERMANOVA: Jaccard Culture-free ----
# isolate endophyte community data 
data.free.comm <- select(data.free, starts_with('Otu'))

# Remove columns with 0 or singletons
data.free.comm[colSums(data.free.comm) > 1] -> data.free.comm

# Hellinger transformation
comm.hell <- decostand(data.free.comm, method = 'total')

# Jaccard dissimilarity
jacc.dist <- vegdist(comm.hell, method = 'jaccard', upper = F,
                     binary = T) 

# Morisita-Horn dissimilarity
horn.dist <- vegdist(comm.hell, method = 'horn', upper = F)

# check for homogeneity of group dispersions
jacc.betadisper <- betadisper(jacc.dist, group = data.all$Genus)
permutest(jacc.betadisper, permutations = 99)
anova(jacc.betadisper)
plot(jacc.betadisper)

jacc.adonis <- adonis2(jacc.dist ~ Native.status * Genus + PCNM1,
                       by = 'term',
                       data = data.free,
                       strata = data.free$site)
jacc.adonis

# write out results
as.data.frame(jacc.adonis) -> jacc.adonis
write.csv(jacc.adonis, 'results/SupplementaryTableS3_CF_Jaccard.csv', row.names = T)

#### NMDS plot: Jaccard Culture-free ----
jacc.mds <- metaMDS(jacc.dist, dist = 'bray',
                    try = 1000, trymax = 1000,
                    tidy = T, shrink = T, k = 2)
jacc.stress <- jacc.mds$stress

#### Make data frame for plot ----
data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                          NMDS2 = jacc.mds$points[,2],
                          sample.no = data.free$SampleID,
                          # site = data.free$site,
                          Native.status = data.free$Native.status,
                          Genus = data.free$Genus, 
                          Species = data.free$Species)

jacc.plot.s4a <- ggplot(data = data.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(aes(fill = Genus),
                 shape = 21,
                 color = 'black',
                 size = 2) +  # Using 'Species' for shape too, or adjust as needed
  stat_ellipse(aes(group = Native.status,
                   color = Native.status),  # map color to Native.status
               geom = "polygon",
               linetype = "dashed",
               fill = NA,
               linewidth = 0.8) +
  scale_fill_brewer(type = 'div') +
  scale_color_manual(values = c('darkgrey', 'black')) +    
  # geom_text(aes(label = sample.no),  # Change to the desired label column
  #           vjust = -1,              # Adjust vertical position
  #           size = 3) +              # Adjust text size
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = 'Black'),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = 'right')

# Supplementary Figure S4 ----

jacc.plot.s4a +
  jacc.plot.s4b +
    plot_layout(
    nrow = 1,
    ncol = 2,
    guides = "collect"
  ) + 
  plot_annotation(tag_levels = "A") -> SupplFig.S4

ggsave("figure/SupplementaryFigureS4.jpeg",
       plot = SupplFig.S4,
       width = 10, height = 6,
       device = 'jpeg',
       dpi = 600)

    
# PERMANOVA: Morisita Horn Culture-free ----
# check for homogeneity of group dispersions
horn.betadisper <- betadisper(horn.dist, group = data.free$Genus)
permutest(horn.betadisper, permutations = 99)
anova(horn.betadisper)
plot(horn.betadisper)

horn.adonis <- adonis2(horn.dist ~ Native.status * Genus + PCNM1,
                       by = 'term',
                       data = data.free,
                       strata = data.free$site)
horn.adonis

# write out results
as.data.frame(horn.adonis) -> horn.adonis
write.csv(horn.adonis,
          'results/SupplementaryTableS3_CF_MorisitaHorn.csv',
          row.names = T)

#### NMDS plot: Morisa Horn Culture-free ----
horn.mds <- metaMDS(horn.dist, dist = 'bray',
                    try = 1000, trymax = 1000,
                    tidy = T, shrink = T, k = 2)
horn.stress <- horn.mds$stress

#### Make data frame for plot ----
data.scores <- data.frame(NMDS1 = horn.mds$points[,1],
                          NMDS2 = horn.mds$points[,2],
                          sample.no = data.free$SampleID,
                          Native.status = data.free$Native.status,
                          Genus = data.free$Genus, 
                          Species = data.free$Species)

data.scores$Species <- as.factor(data.scores$Species)

horn.plot.1b <- ggplot(data = data.scores,
       aes(x = NMDS1,
           y = NMDS2)) +
  geom_point(aes(fill = Genus),
             color = "Black",
             shape = 21,
             size = 2) +  # Using 'Species' for shape too, or adjust as needed
  stat_ellipse(aes(group = Native.status,
                   color = Native.status),  # map color to Native.status
               geom = "polygon",
               linetype = "dashed",
               fill = NA,
               linewidth = 0.8) +
  scale_fill_brewer(type = 'div') +
  scale_color_manual(values = c('darkgrey', 'black')) +    
  # geom_text(aes(label = sample.no),  # Change to the desired label column
  #           vjust = -1,              # Adjust vertical position
  #           size = 3) +              # Adjust text size
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = 'Black'),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = 'right')

# Figure 1 -----
# For creation of Figure 1. Have to run analysis.Q1_Metabolome.R prior to 
# running the beow code as pca1 comes from there.
pca1 +
  horn.plot.1b +
  horn.plot.1c +
    plot_layout(
    nrow = 1,
    ncol = 3,
    guides = "collect"
  ) + 
  plot_annotation(tag_levels = "A") -> Fig.1

ggsave("figure/Figure1.jpeg",
       plot = Fig.1,
       width = 16, height = 5,
       device = 'jpeg',
       dpi = 600)

# Taxonomic analyses ----
data.free.tax <- read.csv('data/G045.TaxonomicData_CultureFree.csv')
data.cult.tax <- read.csv('data/G045.TaxonomicData_CultureBased.csv')

### Culture-free ----
# Reformat data for taxonomic analyses
data.free %>%
  select(SampleID, Genus, Native.status, starts_with('Otu')) %>%
  gather(key = 'Otu', value = 'Read.Count',
         -SampleID, -Genus, -Native.status) %>%
  filter(Read.Count >= 1) %>%
  rename(Host.genus = Genus) %>%
  mutate(count = 1) -> data.cf.long

# Add taxonomic data
data.free.tax %>%
  right_join(data.cf.long, by = 'Otu') -> data.cf.long

# Remove low abundance groups (less than 10 OTU and less than 100 reads)
data.cf.long %>%
  filter(Phylum %in% c('Ascomycota', 'Basidiomycota',
                       'Unknown', 'Zoopagomycota')) -> data.cf.filtered

### Culture-based ----
# Reformat data for taxonomic analyses
data.cult.short %>%
  select(SampleID, Genus, Native.status, starts_with('cOtu')) %>%
  gather(key = 'cOtu', value = 'Occurrence',
         -SampleID, -Genus, -Native.status) %>%
  mutate(count = 1) %>%
  rename(Host.genus = Genus) %>%
  filter(Occurrence >= 1) -> data.cb.long

# Add taxonomic data
data.cult.tax %>%
  rename(cOtu = 'OTU.tbas') %>%
  distinct(cOtu, Phylum, Class, Order) %>%
  right_join(data.cb.long, by = 'cOtu') -> data.cb.long

# Remove low abundance groups
data.cb.long %>%
  filter(Phylum %in% c('Ascomycota', 'Basidiomycota')) -> data.cb.filtered

## Assessment of taxonomy as a function of Genus: Culture-free ----

#### Phylum ----
data.cf.filtered %>%
  group_by(Host.genus, Phylum) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Phylum, value = Total.count, fill = 0) -> data.tax.phylum
data.tax.phylum <- as.data.frame(data.tax.phylum)
row.names(data.tax.phylum) <- data.tax.phylum$Host.genus
data.tax.phylum.short <- data.tax.phylum[-1]

chisq.test(data.tax.phylum.short) # Significant

# Remove temporary data files
rm(data.tax.phylum); rm(data.tax.phylum.short)

#### Class ----
data.cf.filtered %>%
  group_by(Host.genus, Class) %>%
  summarise(Total.count = sum(count)) %>%
  filter(Class %in% c('Agaricomycetes', 'Dothideomycetes', 'Pezizomycetes',
                      'Saccharomycetes', 'Sordariomycetes',
                      'Unknown', 'Ustilaginomycetes', 'Zoopagomycetes')) %>%
  spread(key = Class, value = Total.count, fill = 0) -> data.tax.class

# Replace NAs with 0
data.tax.class <- as.data.frame(data.tax.class)
rownames(data.tax.class) <- data.tax.class$Host.genus
data.tax.class.short <- data.tax.class[-1]

data.tax.class.short <- select(data.tax.class.short,
                               where(~ sum(.x, na.rm = TRUE) >= 3))

chisq.test(data.tax.class.short) # Significant

rm(data.tax.class); rm(data.tax.class.short)

# Plot

data.cf.filtered %>%
  group_by(Host.genus, Class) %>%
  filter(Class != 'Zoopagomycetes') %>%
  summarise(Total.count = sum(count)) %>%
  filter(Total.count >= 3) -> data.cf.phyla

# Group Classes by phyla
data.cf.phyla$Class <- factor(data.cf.phyla$Class,
       levels = c('Agaricomycetes', 'Tremellomycetes', 'Ustilaginomycetes',
                  'Dothideomycetes', 'Eurotiomycetes', 'Pezizomycetes',
                  'Saccharomycetes', 'Sordariomycetes','Unknown'))

Supplfig.S5a <- ggplot(data.cf.phyla,
       aes(x = Host.genus,
             y = Total.count,
             fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "OTU occurrence") +
  scale_fill_manual(values = c('#a50026','#fdae61','#ffffbf',
                               '#e0f3f8','#abd9e9','#74add1','#4575b4','#313695',
                               'darkgrey')) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_blank())

### Ascomycota orders ----
data.cf.long %>%
  filter(Phylum == 'Ascomycota') %>%
  group_by(Host.genus, Order) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Order, value = Total.count, fill = 0) -> data.tax.asco.order

# Organize dataframe
data.tax.asco.order <- as.data.frame(data.tax.asco.order)
rownames(data.tax.asco.order) <- data.tax.asco.order$Host.genus
data.tax.asco.order.short <- data.tax.asco.order[-1]

data.tax.asco.order.short <- select(data.tax.asco.order.short,
                               where(~ sum(.x, na.rm = TRUE) >= 3))

chisq.test(data.tax.asco.order.short) # No significant difference

rm(data.tax.asco.order); rm(data.tax.asco.order.short)

# Plot
Supplfig.S5b <- data.cf.filtered %>%
  group_by(Host.genus, Order) %>%
  filter(Phylum == 'Ascomycota') %>%
  summarise(Total.count = sum(count)) %>%
  filter(Total.count >= 2) %>%
  ggplot(aes(x = Host.genus,
             y = Total.count,
             fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('lightgrey','#f7fbff','#deebf7','#c6dbef','#9ecae1',
                               '#6baed6','#4292c6','#2171b5','#08519c','#08306b',
                               'darkgrey')) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_blank())

### Basidiomycota order----
data.cf.long %>%
  filter(Phylum == 'Basidiomycota') -> data.Basid

data.Basid %>%
  group_by(Host.genus, Order) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Order, value = Total.count, fill = 0) -> data.Basid.order

# Organize dataframe
data.Basid.order <- as.data.frame(data.Basid.order)
rownames(data.Basid.order) <- data.Basid.order$Host.genus
data.Basid.order.short <- data.Basid.order[-1]

data.Basid.order.short <- select(data.Basid.order.short,
                               where(~ sum(.x, na.rm = TRUE) >= 3))

chisq.test(data.Basid.order.short) # Significant difference

rm(data.Basid.order, data.Basid.order.short, data.Basid)

# Plot
Supplfig.S5c <- data.cf.filtered %>%
  group_by(Host.genus, Order) %>%
  filter(Phylum == 'Basidiomycota',
         Class %in% c('Agaricomycetes', 'Tremellomycetes', 'Ustilaginomycetes')) %>%
  summarise(Total.count = sum(count)) %>%
  filter(Total.count >= 3) %>%
  ggplot(aes(x = Host.genus,
             y = Total.count,
             fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "") +
  scale_fill_manual(values = c('#990000','#d7301f','#ef6548','#fc8d59','#fdbb84',
                    '#fdd49e','#fee8c8','#8c2d04')) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_blank())

## Assessment of taxonomy as a function of Genus: Culture-based ----

#### Phylum ----
data.cb.filtered %>%
  group_by(Host.genus, Phylum) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Phylum, value = Total.count, fill = 0) -> data.tax.phylum
data.tax.phylum.short <- data.tax.phylum[-1]
data.tax.phylum.short <- select(data.tax.phylum.short,
                               where(~ sum(.x, na.rm = TRUE) >= 1))

chisq.test(data.tax.phylum.short) # Not significant

# Remove temporary data files
rm(data.tax.phylum); rm(data.tax.phylum.short)

#### Class ----
data.cb.filtered %>%
  group_by(Host.genus, Class) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Class, value = Total.count, fill = 0) -> data.tax.class

# Replace NAs with 0
data.tax.class.short <- data.tax.class[-1]

chisq.test(data.tax.class.short) # Not significant difference

rm(data.tax.class); rm(data.tax.class.short)

# Plot
data.cb.filtered %>%
  group_by(Host.genus, Class) %>%
  summarise(Total.count = sum(count)) -> data.cb.phylum

# Group Classes by phyla
data.cb.phylum$Class <- factor(data.cb.phylum$Class,
       levels = c('Agaricomycetes','Cystobasidiomycetes','Microbotryomycetes','Tremellomycetes',
                  'Ustilaginomycetes',
                  'Dothideomycetes', 'Eurotiomycetes', 'Lichinomycetes',
                  'Sordariomycetes'))

# Pezizo: #74add1, Saccharo #4575b4
Supplfig.S5d <- ggplot(data.cb.phylum, 
       aes(x = Host.genus,
             y = Total.count,
             fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "OTU occurrence") +
  scale_fill_manual(values = c('#a50026','#d73027', '#f46d43','#fdae61','#ffffbf',
                               '#e0f3f8','#abd9e9','#74add1','#313695')) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_blank())

### Ascomycota orders ----
data.cb.filtered %>%
  filter(Phylum == 'Ascomycota') %>%
  group_by(Host.genus, Order) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Order, value = Total.count, fill = 0) -> data.tax.asco.order

data.tax.asco.order.short <- data.tax.asco.order[-1]

chisq.test(data.tax.asco.order.short) # Significant difference

rm(data.tax.asco.order); rm(data.tax.asco.order.short)

# Plot - Order
Supplfig.S5e <- data.cb.filtered %>%
  filter(Phylum == 'Ascomycota') %>%
  group_by(Host.genus, Order) %>%
  summarise(Total.count = sum(count)) %>%
  filter(Total.count >= 2) %>%
  ggplot(aes(x = Host.genus,
             y = Total.count,
             fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "") +
  scale_fill_manual(values = c('#deebf7','#c6dbef','#9ecae1',
                               '#6baed6','#4292c6','#2171b5','#08519c','#08306b')) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_blank())

### Basidiomycota order----
data.cb.filtered %>%
  filter(Phylum == 'Basidiomycota') -> data.Basid

data.Basid %>%
  group_by(Host.genus, Order) %>%
  summarise(Total.count = sum(count)) %>%
  spread(key = Order, value = Total.count, fill = 0) -> data.Basid.order

data.Basid.order.short <- data.Basid.order[-1]
data.Basid.order.short <- select(data.Basid.order.short,
                               where(~ sum(.x, na.rm = TRUE) >= 1)
                               )

chisq.test(data.Basid.order.short) # No significant difference

rm(data.Basid.order, data.Basid.order.short, data.Basid)

# Plot - Order
Supplfig.S5f <- data.cb.filtered %>%
  filter(Phylum == 'Basidiomycota') %>%
  group_by(Host.genus, Order) %>%
  summarise(Total.count = sum(count)) %>%
  # filter(Total.count >= 2) %>%
  ggplot(aes(x = Host.genus,
             y = Total.count,
             fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "") +
  scale_fill_manual(values = c('#990000','#d7301f','#ef6548','#fc8d59','#fdbb84',
                    '#fdd49e','#fee8c8','#8c2d04')) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = 'black',
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_blank())

# Supplementary Figure S5 ----
Supplfig.S5a +
  Supplfig.S5b +
  Supplfig.S5c +
  Supplfig.S5d +
  Supplfig.S5e +
  Supplfig.S5c +
    plot_layout(
    nrow = 2,
    ncol = 3
  ) + 
  plot_annotation(tag_levels = "A") -> SupplFig.S5

ggsave("figure/SupplementaryFigureS5.jpeg",
       plot = SupplFig.S5,
       width = 14, height = 10,
       device = 'jpeg',
       dpi = 600)
