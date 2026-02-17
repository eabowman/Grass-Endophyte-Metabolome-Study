##################################################################################
##################################################################################
# ------------------------------------------------------------------------------ #
# ////// Question 3 Script. LASSO bioactive metabolites - OTU interactions ///// #
# ------------------------------------------------------------------------------ #
##################################################################################
##################################################################################
# Colin Richard Morrison
# The University of Texas at Austin, Department of Integrative Biology 
# Brackenridge Field Laboratory, 2907 Lake Austin Blvd., Austin TX 78703
# *crmorrison@utexas.edu
# additional analyses for questions 1 and 3 paper by Dr. Liz Bowman and Colin Morrison accompany this script on GitHub

# Goal: Identify  metabolite–fungus associations to guide hypotheses about chemical mediation on host plants


setwd("")

library(ggplot2)
library(vegan)
library(tidyverse)
library(dplyr)
library(reshape2)
library(glmnet)
library(bipartite)


##############################################################################
### --- 1. Wrangle and format final OTU-host matrix for LASSO analysis --- ###
##############################################################################
### load the data 
fe<-read.csv("G045.CF.OTU.long.csv")
dim(fe) # 230 OTUs

### Identify OTUs that are present >1 SampleID 
# and with > 6 reads on all CF samples, 
# or all reads abundancesfor CB samples
fe_filtered_no_singles <- fe %>%
  filter(occurence >= 6) %>%                # read abundance threshold
  group_by(Otu,Phylum,Class,Order) %>%
  filter(n_distinct(SampleID) > 1) %>%      # keep only OTUs present in >1 sample
  ungroup()

### Visualize the number of OTUs with different read abundances 
# Summarize number of unique hosts per OTU 
otu_host_counts <- fe_filtered_no_singles %>%
  group_by(Otu) %>%
  summarise(n_hosts = n_distinct(SampleID)) %>%
  ungroup()

# plot histogram
hist<-ggplot(otu_host_counts, aes(x = n_hosts)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue", boundary = 0.5) +
  scale_x_continuous(breaks = seq(1, max(otu_host_counts$n_hosts), 1),
                     expand = expansion(mult = c(0.15, 0.05))) +
  labs(
    x = "Number of individual grass samples",
    y = "Number of OTUs",
    title = "Culture-free fungal endophyte host specificity") +
  theme(
    panel.grid.major.x = element_line(colour = "grey85", size=0.25),
    panel.grid.major.y = element_line(colour = "grey85", size=0.25),
    panel.background   = element_rect(size = 0.8, linetype = 'solid',
                                      colour = "black", fill="white"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust=0,vjust=2.5,face="bold",  size=22,
                              margin = margin(t = 5)),
    axis.title.x = element_text(margin = margin(t = 30),
                                vjust=1.5, size=16),
    axis.title.y = element_text(margin = margin(t = 30),
                                vjust=9.0, size=16),
    axis.text.x  = element_text(hjust=0.5, size=14, color="black"),
    axis.text.y  = element_text(size=14, color="black"),
    plot.margin  = margin(t = 10, r = 40, b = 25, l = 30, unit = "pt")
  )
hist # 
  
### 'cast' into a molten food web network with rows of metabolites and columns of FE OTUs
wide_lasso <-dcast(
  fe_filtered_no_singles, # 
  SampleID ~ Otu,
  value.var="occurence",
  fun.aggregate=sum
) 
# clean up the OTU-host network
row.names(wide_lasso)<-wide_lasso$SampleID 
wide_lasso <-wide_lasso[,-c(1)]
OTUs = t(wide_lasso)  # rows are  FE OTUs, columns are grass samples, entries are OTU relative abundances
OTUs <- as.data.frame(OTUs)


### --- set .Rdata object for saving LASSO  outputs
date = 20251111 # date analysis run
class = "total" # specify which metabolomic sube=set
fungi = "CF_>=5reads_>=1grass_samples" # specify OTU abundance thresholh
outfile = paste("metabolomics/results/LASSO_v3_",class,"_",fungi,"_",date, ".RData",sep="")


###########################################
### --- 2. Perform LASSO regression --- ###
###########################################
### make metabolite classification dataframe for LASSO results below
metab.master<-read.csv("G045_v2_metabolmes_total_20251018.csv")
head(metab.master)
# subset chemotaxonomic data for LASSO below
siri_real_lasso=metab.master[,c("X","row.ID","row.retention.time","row.m.z",
                                "NPC.pathway","NPC.superclass","NPC.class","ClassyFire.most.specific.class")]
rownames(siri_real_lasso)=siri_real_lasso$row.id
siri_real_lasso=siri_real_lasso[,-c(1)]


###  make metabolite abundance dataframe for LASSO
metab<-load("metabolomics/G045_v2_jaccard_total_classified_90_pathway_prob_20251017.RData")
dim(heat.class) #1723 compounds that were classified to biochemical pathway >= 90% confidence by SIRIUS

# rows must be metabolites, columns must be host plant samples
master_real=heat.class ### make sure to order the host names same as herbivore dataset
# reorder into ABC order 
master_real_lasso <- master_real[ , order(names(master_real))] 
# subset only grass samples with metabolite abundances
names(master_real_lasso)
master_real_lasso=master_real_lasso[,c(1:36)]


### --- make new objects for LASSO for-loop
metab.plants = master_real_lasso 
metab_cols <- colnames(metab.plants)

otu_cols <- colnames(OTUs)

# remove SampleID column(s) not present in both data trables
missing_in_otu <- setdiff(metab_cols, otu_cols)

metab.plants <- metab.plants[, !(colnames(metab.plants) %in% missing_in_otu)]
metab.plants <- metab.plants[, otu_cols]
identical(colnames(OTUs), colnames(metab.plants)) # [1] TRUE


# make empty dataframe for Lasso loop to fill
coef.metab.OTUs = as.data.frame(matrix(0,nrow = nrow(metab.plants), ncol = nrow(OTUs)))
names(coef.metab.OTUs) = row.names(OTUs)
row.names(coef.metab.OTUs) = row.names(metab.plants) # rows are compounds, columns are endophytes

# count how many OTUs there are and update 
dim(OTUs) # 44 35 (>5 reads & >1 grasssample)

# define number of OTUs that will be abalyed
nendo = 44 

### for-loop to calculate regression coefficients associating grass metabolites with OTU presence begins here
reclist = rep(0, length = nendo)

for(k in 1:nendo){
  # Check if the index is valid before proceeding
  if (k > nrow(OTUs)) {
    cat("Index k =", k, "is out of bounds for the OTUs matrix.\n")
    next
  }
  # Try to extract and check the data for the current OTU
  tryCatch({
    y_k = t(OTUs[k,])
    # Check for both variation and non-zero entries
    if (length(unique(y_k)) > 1 && sd(y_k) > 0) {
      cvfit = cv.glmnet(x = t(metab.plants), y = y_k) # Poisson dist is default distribution
      coef_obj = coef(cvfit, s = "lambda.min")
      compeffs = rownames(coef_obj)[which(coef_obj[,1] != 0)]
      
      if(length(compeffs) > 1){
        cat("Showing compounds of effect for otu ", k, "\n")
        reclist[k] = 1
        
        non_zero_coeffs = coef_obj[which(coef_obj[,1] != 0), ]
        non_zero_coeffs = non_zero_coeffs[which(names(non_zero_coeffs) != "(Intercept)")]
        
        for(j in 1:length(non_zero_coeffs)){
          row_name = names(non_zero_coeffs)[j]
          if(row_name %in% row.names(coef.metab.OTUs)){
            coef.metab.OTUs[row_name, k] = non_zero_coeffs[j]
          }
        }
      }
    } else {
      cat("Skipping otu ", k, " because its response vector is constant or has insufficient variation (y_k values are:", toString(y_k), ").\n")
    }
  }, error = function(e) {
    # This block will be executed if an error occurs
    cat("An error occurred for OTU", k, ":", conditionMessage(e), "\n")
    cat("Problematic y_k data for k =", k, ":\n")
    print(y_k)
  })
}

# number of OTUs analyzed:
k
# number of OTUs with significant regression coefficients (= associated w bioactive metabolites):
sum(reclist)
# ratio of OTUs with bioactive metabolites:
cat(sum(reclist), " out of ", k, " OTUs had a compound of effect", "\n")
# confirm how many bioactive compounds:
length(which(rowSums(coef.metab.OTUs) != 0)) 



### --- calculate net value of bioactive  metabolites 
dim(coef.metab.OTUs) 

nonzero_metabs <- coef.metab.OTUs[rowSums(coef.metab.OTUs) != 0, ]
length(which(rowSums(nonzero_metabs) < 0)) # negatively associated metabs
length(which(rowSums(nonzero_metabs) > 0)) # positively associated metabs

### --- mean of nonzero coefficients only, plus inversion so strongest repellent = +1 and strongest attractant = –1
HAM.coef.metab.OTUs <- data.frame(
  metabolite = rownames(nonzero_metabs),
  mean_coef = apply(nonzero_metabs, 1, function(x) mean(x[x != 0]))
)

### --- log10 transform magnitudes 
HAM.coef.metab.OTUs$log_coef <- log10(abs(HAM.coef.metab.OTUs$mean_coef) + 1e-10)

### --- rescale magnitudes between 1 and 2  
HAM.coef.metab.OTUs$scaled_coef <- HAM.coef.metab.OTUs$log_coef / max(HAM.coef.metab.OTUs$log_coef)

### --- restore original sign to preserve biological direction of the effect on OTUs
HAM.coef.metab.OTUs$scaled_coef <- sign(HAM.coef.metab.OTUs$mean_coef) * HAM.coef.metab.OTUs$scaled_coef

### --- add NPclassifier annotations to the summary
# Convert rownames of siri_real_lasso into a column to merge on
siri_real_lasso$metabolite <- rownames(siri_real_lasso)
siri_annot <- siri_real_lasso[, c("metabolite","NPC.pathway","NPC.superclass",
                                  "NPC.class", "ClassyFire.most.specific.class")]
HAM.coef.metab.OTUs.annot <- merge(HAM.coef.metab.OTUs,siri_annot,by = "metabolite",all.x = TRUE)

### --- save lasso data
save(coef.metab.OTUs, OTUs, fe, wide_lasso,
     heat.class,master_real, master_real_lasso, metab.plants,
     metab.master, siri_real_lasso, nonzero_metabs,
     HAM.coef.metab.OTUs, HAM.coef.metab.OTUs.annot, 
     siri.pos, siri.neg,  file = outfile)




#####################################################################################
### --- 3. create bipartite food webs OTUs interacting with LASSO metabolites --- ###
#####################################################################################
### load the data
lasso_data<-load("metabolomics/results/LASSO_v3_total_CF_>=6reads_>=2grass_samples_20251105.RData")
# FINAL Culture-free data: "LASSO_v3_total_CB_>=1reads_>=2grass_samples_20251111.RData"
# Final culture-based data: "LASSO_v3_total_CF_>=6reads_>=2grass_samples_20251105.RData"


###  Add NPC.pathway to LASSO metabolites data
annot <- HAM.coef.metab.OTUs.annot %>%
  dplyr::select(metabolite, NPC.pathway, NPC.superclass, ClassyFire.most.specific.class) %>%
  mutate(
    metabolite = as.character(metabolite),
    NPC.pathway = as.character(NPC.pathway)
  ) %>%
  filter(!is.na(NPC.pathway) & NPC.pathway != "")


# put LASSO metabolite_OTU pairs into long format
lasso_links <- nonzero_metabs %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Metabolite") %>%
  tidyr::pivot_longer(
    cols = -Metabolite,
    names_to = "OTUs",
    values_to = "coef"
  ) %>%
  dplyr::mutate(
    OTUs = as.character(OTUs),
    Metabolite = as.character(Metabolite)
  ) %>%
  dplyr::filter(!is.na(coef) & coef != 0) %>%
  dplyr::select(Metabolite, OTUs)


### put OTUs abundance data in long format
OTUs_long <- OTUs %>%                
  as.data.frame() %>%
  tibble::rownames_to_column(var = "OTUs") %>%
  pivot_longer(
    cols = -OTUs,
    names_to = "SampleID",
    values_to = "OTUs_abund"
  ) %>%
  mutate(
    OTUs = as.character(OTUs),
    SampleID = as.character(SampleID),
    OTUs_abund = as.numeric(OTUs_abund)
    )

### filter to only the grass samples where OTUs present 
OTUs_present <- OTUs_long %>%
  filter(OTUs_abund > 0)

### join LASSO pairs to OTU presence by OTUs
interactions <- lasso_links %>%
  inner_join(OTUs_present, by = "OTUs")   # join column is "OTUs" in both tables

### attach NPC.pathway annotations & drop pathways without annotation (= NA)
interactions_annot <- interactions %>%
  mutate(Metabolite = as.character(Metabolite)) %>%
  left_join(annot, by = c("Metabolite" = "metabolite")) %>%
  filter(!is.na(NPC.pathway) & NPC.pathway != "")


### --- Summarize associations by NPC.pathway × OTUs and annotate the OTU taxonomy
fe<-read.csv("G045.CF.OTU.long.csv")
head(fe)
dim(fe) # 

### Identify OTUs that are present >1 SampleID 
# and with > 6 reads on all CF samples, 
# or all reads abundancesfor CB samples
fe_filtered_no_singles <- fe %>%
  filter(occurence >= 6) %>%                # first filter for >=11 occurrences
  group_by(Otu,Phylum,Class,Order) %>%
  filter(n_distinct(SampleID) > 1) %>%      # keep only OTUs present in >=2 samples
  ungroup()


### calculate the number unique metabolites from each NPC.pwathway associated with OTUs
network_summary <- interactions_annot %>%
  group_by(NPC.pathway, OTUs) %>%
  summarise(
    n_metabolites = n_distinct(Metabolite), # Number unique  metabolites in a pathway selected by LASSO and were co-present with that OTU >=2 grass samples
    n_samples_with_OTU = n(), # Number grass samples support the interaction between a given NPC.pathway and a specific OTU
    OTU_interactions = n_distinct(OTUs),
    .groups = "drop"
  )

### add OTU taxonomy to the dataframe
network_summary_tax <- network_summary %>%
  left_join(
    fe_filtered_no_singles %>%
      distinct(Otu, Phylum, Class, Order),
    by = c("OTUs" = "Otu")
  ) %>%
  group_by(Phylum, Class, Order) %>%
  mutate(
    Phylum_unique = paste0(Phylum, "_", row_number()),
    Class_unique = paste0(Class, "_", row_number()),
    Order_unique = paste0(Order, "_", row_number())
  ) %>%
  ungroup()

### remove Choanozoan flagellate from CF dataset (not an FE)
# network_summary_tax <- network_summary_tax[network_summary_tax$OTUs != "Otu413", ]


### final number of metabolites that are included in the network:
dim(network_summary) # [1] 27


### --- Cast to wide form NPC.pathway × OTUs food-web like matrix
wide <- dcast(
  network_summary_tax,
  NPC.pathway ~ Phylum_unique,
  value.var = "n_samples_with_OTU", #  
  fun.aggregate = sum,
  fill = 0
)
rownames(wide) <- wide$NPC.pathway
wide <- as.matrix(wide[ , -1, drop = FALSE])


### customize the order of OTUs
# CB FEs (>= 1 reads)
new_col_order <- c("Ascomycota_1","Ascomycota_2","Ascomycota_3","Ascomycota_4","Ascomycota_5",
                   "Ascomycota_6","Ascomycota_7","Ascomycota_8","Ascomycota_9","Ascomycota_10",
                   "Ascomycota_11","Ascomycota_12","Ascomycota_13","Ascomycota_14",
                   "Basidiomycota_1","Basidiomycota_2")
### CF FEs (>= 6 reads)
new_col_order <- c(
  "Ascomycota_1", "Ascomycota_2", "Ascomycota_3", "Ascomycota_4", "Ascomycota_5",
  "Basidiomycota_1", "Basidiomycota_10", "Basidiomycota_11", "Basidiomycota_12", "Basidiomycota_2",
  "Basidiomycota_3", "Basidiomycota_4", "Basidiomycota_5", "Basidiomycota_6", "Basidiomycota_7",
  "Basidiomycota_8", "Basidiomycota_9", "Unknown_1", "Unknown_2",
  "Unknown_3", "Unknown_4", "Unknown_5", "Unknown_6", "Unknown_7")

wide <- wide[, new_col_order, drop = FALSE] 




### --- plot it
#pdf("metabolomics/plots/bipartite_lasso_hams_CB_>=1reads_>2grass_samples.pdf", width = 11, height = 9)

plotweb(
  wide,
  method = "normal",
  text.rot = 90,
  low.y = 0.5,    
  high.y = 1.4,      
  col.interaction = "grey50",  
  col.high = "grey20",                 
  col.low = "grey20",
  y.width.low = 0.02, 
  y.width.high = 0.02,
  y.lim = c(0.03, 1.68),
  empty = TRUE)
mtext(text="B ",adj = 0.01,
      side = 3, line = -1.5, outer = FALSE,cex=1.5)
mtext(text="Culture-based fungal endophytes",
      side = 3, line = -1.5, outer = FALSE,cex=1.7)
mtext(text="Grass biochemical pathways", side = 1, line = -1.0, outer = FALSE,cex=1.7)
mtext(text="(N = 35 bioactive metabolites)", side = 1, line = 0.2, outer = FALSE,cex=1.2)

#dev.off()




###########################################################################
### --- 4. make OTU-LASSO metabolite NPC.pathway and NPC.superclass --- ###
###########################################################################

### NPC biochemical pathway of origin
pathway_summary<-interactions_annot %>%
  filter(!is.na(NPC.pathway) & NPC.pathway != "") %>%
  group_by(NPC.pathway) %>%
  summarise(
    n_metabolites = n_distinct(Metabolite),
    OTU_interactions = n_distinct(OTUs)  
    ) %>%
  mutate(
    percent_total = 100 * n_metabolites / sum(n_metabolites),
    percent_total = sprintf("%.2f", percent_total)
  ) %>%
  arrange(desc(as.numeric(percent_total))) %>%
  dplyr::select(
    NPC.pathway,
    n_metabolites,
    percent_total,
    OTU_interactions)

pathway_summary


### NPC.superclass
# Define desired pathway order (top → bottom in plotweb)
pathway_order <- c(
  "Shikimates and Phenylpropanoids",
  "Terpenoids",
  "Fatty acids",
  "Carbohydrates",
  "Polyketides",
  "Amino acids and Peptides",
  "Alkaloids")

superclass_summary <- interactions_annot %>%
  filter(!is.na(NPC.superclass) & NPC.superclass != "") %>%
  group_by(NPC.superclass) %>%
  summarise(
    n_metabolites = n_distinct(Metabolite),
    OTU_interactions = n_distinct(OTUs),
    NPC.pathway = paste(unique(na.omit(NPC.pathway)), collapse = "; "),
    .groups = "drop",
    CF.most.specific = paste(unique(na.omit(ClassyFire.most.specific.class)), collapse = "; ")
  ) %>%
  mutate(
    percent_total_num = 100 * n_metabolites / sum(n_metabolites), # numeric version
    percent_total = sprintf("%.2f", percent_total_num),           # formatted version
    NPC.pathway = factor(NPC.pathway, levels = pathway_order)     # enforce pathway ordering
  ) %>%
  arrange(NPC.pathway, desc(percent_total_num)) %>%  # <-- secondary ordering
  dplyr::select(
    NPC.pathway,
    NPC.superclass,
    CF.most.specific,
    n_metabolites,
    percent_total,
    OTU_interactions
  )
superclass_summary


### ---- save lasso data once more
save(coef.metab.OTUs, OTUs, fe, wide_lasso,
     heat.class,master_real, master_real_lasso, metab.plants,
     metab.master, siri_real_lasso, superclass_summary, pathway_summary,
     HAM.coef.metab.OTUs, HAM.coef.metab.OTUs.annot, 
     siri.pos, siri.neg,  file = outfile)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################