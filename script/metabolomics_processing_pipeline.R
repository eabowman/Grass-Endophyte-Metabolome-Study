###############################################################
###############################################################
# ----------------------------------------------------------- #
# ////// Grass Metabolomics Data Processing (G045) ///// #
# ----------------------------------------------------------- #
###############################################################
###############################################################

# Colin Richard Morrison*
# The University of Texas at Austin, Department of Integrative Biology 
# Brackenridge Field Laboratory, 2907 Lake Austin Blvd., Austin TX 78703
# *crmorrison@utexas.edu


library(phytools)
library(picante)
library(caper)
library(vegan)
library(MASS)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(xlsx)


getwd()
setwd("")


### --- samples metadata table
meta = read.csv("G045_metadata_editedEAB.csv", header = T)
head(meta)
nrow(meta)

### --- MZmine features tables (ion counts in the samples)
feat = read.csv("G045_features_quantification_table", sep = "") 
names(feat)
feat=feat[,-41]  # remove blank columns that remains
head(feat)

### --- Classyfier results
sirid <- read.csv("G045_structure_identifications.csv")
head(sirid)

### --- NPclassifier results
sirca <- read.csv("G045_canopus_formula_summary.csv")
head(sirca)

### --- Network edges GNPS file for chemical similarity (cosine scores)
network = read.delim("G045_20k_20221117.selfloop",
                     header = T, comment.char = "")
head(network)

### --- Save the LC metadata for df assembly later
head(feat)
featLC = feat[,1:3] 

featLC$smiles = NA
featLC$NPC.pathway = NA
featLC$NPC.superclass = NA
featLC$NPC.class = NA
featLC$ClassyFire.most.specific.class = NA

### --- Filter out features with low congidence to be assigned to a biochemical pathway
hist(sirca$NPC.class.Probability)
dim(sirca) # 
sirca = sirca[which(sirca$NPC.pathway.Probability > 0.9),] # do not select this if you want to analyze all detectable metabolites
dim(sirca) #

### --- Populate feat dataframe with LC features, Classyfire, and NPclassifier data
dim(feat) # [1] TX: 8991   46 ---###--- KE: 7810   39  
for(i in 1:nrow(feat)){
  tempid = feat$row.ID[i]
  featLC$smiles[i] = sirid$smiles[which(sirid$mappingFeatureId == tempid)][1]
  featLC$NPC.pathway[i] = sirca$NPC.pathway[which(sirca$mappingFeatureId == tempid)][1]
  featLC$NPC.superclass[i] = sirca$NPC.superclass[which(sirca$mappingFeatureId == tempid)][1]
  featLC$NPC.class[i] = sirca$NPC.class[which(sirca$mappingFeatureId == tempid)][1]
  featLC$ClassyFire.most.specific.class[i] = sirca$ClassyFire.most.specific.class[which(sirca$mappingFeatureId == tempid)][1]
}

dim(feat) 
dim(featLC)
head(featLC)

### --- Make master table with LC features, compound predictions, and ion counts in each sample
mastertable = cbind(featLC,feat[4:ncol(feat)])
dim(mastertable)
# 
head(mastertable)

### clean up samples names 
cols_to_modify <- grepl(".Peak.area", colnames(mastertable))  # Logical vector
colnames(mastertable)[cols_to_modify] <- gsub(".Peak.area", "", colnames(mastertable)[cols_to_modify])
names(mastertable)
### remove 6 digit sample numbers (if present in data)
names(mastertable)
mastertable <- mastertable %>% rename( "144749"= "X01.mzXML",
                                   "144750"= "X02.mzXML",
                                   "144751"= "X03.mzXML",
                                   "144764"= "X04.mzXML",
                                   "144765"= "X05.mzXML",
                                   "144766"= "X06.mzXML",
                                   "144770"= "X07.mzXML",
                                   "144771"= "X08.mzXML",
                                   "144772"= "X09.mzXML",
                                   "144773"= "X10.mzXML",
                                   "144774"= "X11.mzXML",
                                   "144775"= "X12.mzXML",
                                   "144752"= "X13.mzXML",
                                   "144753"= "X14.mzXML",
                                   "144754"= "X15.mzXML",
                                   "144758"= "X16.mzXML",
                                   "144759"= "X17.mzXML",
                                   "144760"= "X18.mzXML",
                                   "144746"= "X19.mzXML",
                                   "144747"= "X20.mzXML",
                                   "144748"= "X21.mzXML",
                                   "144779"= "X22.mzXML",
                                   "144780"= "X23.mzXML",
                                   "144781"= "X24.mzXML",
                                   "144761"= "X25.mzXML",
                                   "144762"= "X26.mzXML",
                                   "144763"= "X27.mzXML",
                                   "144767"= "X28.mzXML",
                                   "144768"= "X29.mzXML",
                                   "144769"= "X30.mzXML",
                                   "144755"= "X31.mzXML",
                                   "144756"= "X32.mzXML",
                                   "144757"= "X33.mzXML",
                                   "144776"= "X34.mzXML",
                                   "144777"= "X35.mzXML",
                                   "144778"= "X36.mzXML")

### --- Select only those compounds that were not found in the blank
dim(mastertable) # 
master.real = mastertable[which(mastertable[,grep("blank.mzXML", names(mastertable))] == 0),]
dim(master.real) #


### now remove the blank column from the df
grep("blank.mzXML", names(mastertable)) # [1] 9 
master.real=master.real[,-c(9)] #
dim(master.real) #




####################################################################################################

### --- make columns identifying total, primary and secondary metabolites 
# shows the different metabolite types in each chemotaxonomic class
levels(as.factor(master.real$ClassyFire.most.specific.class))
levels(as.factor(master.real$NPC.class))
levels(as.factor(master.real$NPC.superclass))
levels(as.factor(master.real$NPC.pathway)) 

### Select a chemical class to analyze in CSCS loop
class = "total_unclassified" ### do not filter by NPC pathway prob on line 68 if you want to analyze all metabolites
class= "total_classified_90_pathway_prob"
class = "secondary_90_pathway_prob"
class = "primary_90_pathway_prob"
# Note: NAs from CLASSY and NPC predictions will not be in metabolite categories 
#       --> So using the following lines only analyzes features that were classifiable by CLASSY and NPC
class = "alkaloids"
class = "peptides"
class = "polyketides"
class = "shikimates"
class = "terpenoids" 
#class = "carbohydrates_90_pathway_prob"
class = "fatty_acids"
class = "tryptophan alkaloids"
class = "flavonoids" 

### --- run this code this before the CSCS loop
heat.itol = master.real # make new object for the upcoming analyses
dim(heat.itol) # (should be same as master.real)

### --- then tun this code for specific chemotaxonomic group you want to analyze
# TOTAL detectable metabolites with LCMS data (this must be same as class object selected above)
if(class == "total_unclassified"){
  heat.class = heat.itol
}

# TOTAL classifiable metabolites through Classyfire or NPClassifier
# make sure to subset compund features >0.90 probability correct classification on L92 first
if(class == "total_classified_90_pathway_prob"){
  heat.class = heat.itol[which(heat.itol$NPC.pathway != "NA"),]
}

# only SECONDARY metabolites 
# make sure to subset compund features >0.90 probability correct classification first
if(class == "secondary_90_pathway_prob"){
  heat.class = heat.itol[which(heat.itol$NPC.pathway %in% c("Alkaloids","Amino acids and Peptides","Polyketides","Shikimates","Shikimates and Phenylpropanoids", "Terpenoids")),]
}

# only PRIMARY metabolites
# make sure to subset compund features >0.90 probability correct classification first
if(class == "primary_90_pathway_prob"){
  heat.class = heat.itol[which(heat.itol$NPC.pathway %in% c("Carbohydrates","Fatty acids")),]
}

# SPECIFIC metabolite classes:
if(class == "alkaloids"){
  heat.class = heat.itol[grep("Alkaloids", heat.itol$NPC.pathway),]
}

if(class == "flavonoids"){
  heat.class = heat.itol[which(heat.itol$NPC.superclass %in% c("Flavonoids","Isoflavonoids")),]
}
dim(heat.class) # [1] 382  43
prob_values <- sirca$NPC.superclass.Probability[match(heat.class$row.ID, sirca$mappingFeatureId)] # Find the corresponding probabilities from sirca
heat.class <- heat.class[which(grepl("flavonoid", heat.class$NPC.superclass, ignore.case = TRUE) & prob_values >= 0.9), ]
dim(heat.class) # [1] 290  43

if(class == "terpenoids"){
  heat.class = heat.itol[grep("Terpenoids", heat.itol$NPC.pathway),]
}

if(class == "tryptophan alkaloids"){
  heat.class = heat.itol[grep("Tryptophan alkaloids", heat.itol$NPC.superclass),]
}


if(class == "peptides"){
  heat.class = heat.itol[grep("Peptides", heat.itol$NPC.pathway),]
}


if(class == "polyketides"){
  heat.class = heat.itol[grep("Polyketides", heat.itol$NPC.pathway),]
}


if(class == "shikimates"){
  heat.class = heat.itol[grep("Shikimates", heat.itol$NPC.pathway),]
}


if(class == "fatty_acids"){
  heat.class = heat.itol[which(heat.itol$NPC.pathway == "Fatty acids"),]
}


#if(class == "carbohydrates"){
#  heat.class = heat.itol[which(heat.itol$NPC.pathway == "Carbohydrates"),]
#}


####################################################################################################
### --- save the master metabolome and classifications as a supplemental table for the paper --- ###
####################################################################################################
write.csv(heat.class,file="G045_metabolmes_total_classified.csv") 



####################################################################################################
####################################################################################################

### --- Remove rows from network data that have compounds that weren't in master.real (blank contaminants) 
### - SOLUTION 1
network.real = network[which(network$CLUSTERID1 %in% heat.class$row.ID),]
nrow(network.real) # [1] 5188 
network.real = network.real[which(network.real$CLUSTERID2 %in% heat.class$row.ID),]
nrow(network.real) # [1] 1443
head(network.real)

### --- Now that the data has been filtered for compounds of interest, remove the LC, and compound prediction data to make it usable for CSCS
names(heat.class)
row.names(heat.class) = heat.class$row.ID
head(heat.class)
heatCSCS = heat.class[,9:44] # 
head(heatCSCS)
names(heatCSCS)

### --- now order the species/samples alphabetically so they line up with phylogenetic distance mtrix later
heatCSCS <- heatCSCS[ , order(names(heatCSCS))]
names(heatCSCS)

sampsByCompounds = as.data.frame(t(as.data.frame(heatCSCS)))
names(sampsByCompounds) = row.names(heatCSCS)
sampsByCompounds[,1:5]
dim(sampsByCompounds)


### --- date that MolecNetsTraits was run
date = 20251023
# if calculating CSCS for species pairs rather than sample pairs, set to TRUE
species = TRUE
outfile = paste("G045_v2_jaccard_", class,"_",date, ".RData",sep="") # change to KE for Kenya metabolomics


## 
nspp = nrow(sampsByCompounds)
ncomps = ncol(sampsByCompounds)
net.comps = c(levels(network.real$CLUSTERID1), levels(network.real$CLUSTERID2))
nspp = nrow(sampsByCompounds)
ncomps = ncol(sampsByCompounds)
pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
names(pairwise.spp) = row.names(sampsByCompounds)
row.names(pairwise.spp) = row.names(sampsByCompounds)
sampsCompsStand = sampsByCompounds
for(i in 1:nrow(sampsByCompounds)){	
  sampsCompsStand[i,] = sampsByCompounds[i,]/sum(sampsByCompounds[i,])
}
diags = pairwise.spp
for (k in 1:nspp){
  sppX = as.character(row.names(sampsCompsStand)[k])
  cat("Comparing ", sppX, " to itself", "\n", sep = "")
  sppXonly = sampsCompsStand[k,which(sampsCompsStand[k,]>0)]
  ncomps = length(sppXonly)
  pairwise.comps = as.data.frame(matrix(0, ncol = ncomps, nrow = ncomps))
  names(pairwise.comps) = names(sppXonly)
  row.names(pairwise.comps) = names(sppXonly)
  for (y in 1:ncomps){
    comp1 = names(sppXonly)[y]
    links.comp1 = network.real[which((network.real$CLUSTERID1 == comp1)|(network.real$CLUSTERID2 == comp1)),]
    for(z in 1:ncomps){
      comp2 = names(sppXonly)[z]
      if(comp2 %in% c(as.character(links.comp1$CLUSTERID1),as.character(links.comp1$CLUSTERID2))){
        if(comp1 == comp2){
          pairwise.comps[y,z] = 1
        }
        if(comp1 != comp2){				
          if((comp1 %in% links.comp1$CLUSTERID1)&(comp2 %in% links.comp1$CLUSTERID2)){
            pairwise.comps[y,z] = network.real$Cosine[which((network.real$CLUSTERID1 == comp1) & (network.real$CLUSTERID2 == comp2))]
          }
          if((comp1 %in% links.comp1$CLUSTERID2)&(comp2 %in% links.comp1$CLUSTERID1)){
            pairwise.comps[y,z] = network.real$Cosine[which((network.real$CLUSTERID2 == comp1) & (network.real$CLUSTERID1 == comp2))]
          }
        }				
      }
    }
  }
  diags[k,k] = sum(((outer(as.numeric(sppXonly), as.numeric(sppXonly)))*pairwise.comps),na.rm = T)	
}
save(sampsByCompounds, pairwise.comps, diags, file = outfile)
for (i in 1:nspp){
  spp1 = as.character(row.names(sampsCompsStand)[i])
  cat("Comparing", spp1, "\n")
  for (j in i:nspp){
    spp2 = as.character(row.names(sampsCompsStand)[j])
    spp1comps = sampsCompsStand[spp1,]
    spp2comps = sampsCompsStand[spp2,]
    spp_pair = rbind(spp1comps,spp2comps)
    paircomps = spp_pair[,which(colSums(spp_pair)>0)]
    ncomps = ncol(paircomps)
    pairwise.comps = as.data.frame(matrix(0, ncol = ncomps, nrow = ncomps))
    names(pairwise.comps) = names(paircomps)
    row.names(pairwise.comps) = names(paircomps)
    for (y in 1:ncomps){
      comp1 = names(paircomps)[y]
      links.comp1 = network.real[which((network.real$CLUSTERID1 == comp1)|(network.real$CLUSTERID2 == comp1)),]
      for(z in 1:ncomps){
        comp2 = names(paircomps)[z]
        if(comp2 %in% c(as.character(links.comp1$CLUSTERID1),as.character(links.comp1$CLUSTERID2))){
          if(comp1 == comp2){
            pairwise.comps[y,z] = 1
          }
          if(comp1 != comp2){				
            if((comp1 %in% links.comp1$CLUSTERID1)&(comp2 %in% links.comp1$CLUSTERID2)){
              pairwise.comps[y,z] = network.real$Cosine[which((network.real$CLUSTERID1 == comp1) & (network.real$CLUSTERID2 == comp2))]
            }
            if((comp1 %in% links.comp1$CLUSTERID2)&(comp2 %in% links.comp1$CLUSTERID1)){
              pairwise.comps[y,z] = network.real$Cosine[which((network.real$CLUSTERID2 == comp1) & (network.real$CLUSTERID1 == comp2))]
            }
          }				
        }
      }
    }
    pairwise.spp[i,j] = pairwise.spp[j,i] = sum(((outer(as.numeric(paircomps[1,]), as.numeric(paircomps[2,])))*pairwise.comps), na.rm = T)/max(diags[i,i], diags[j,j])
  }
}
cscs = pairwise.spp

dim(cscs)

### --- save all the objects 
save(sampsByCompounds,  sampsCompsStand, diags, cscs, heat.class, heat.itol, file = outfile)


### - Bray-Curtis dissimilarity matrix of metabolites (diagonals = 0s)
# compares relative abundances of shared metabolites
braydist = as.data.frame(as.matrix(vegdist(sampsByCompounds, diag = T, upper = T)))
braydist[1:6,1:6]

### - BINARY (presence/absence) Bray-Curtis dissimilarity matrix of metabolites 
braydist.bin = as.data.frame(as.matrix(vegdist(sampsByCompounds, diag = T, upper = T, binary = T)))
dim(braydist.bin) # [1] 89 89
braydist.bin[1:6,1:6]
range(braydist.bin) # 0.0000000 0.9857973
names(braydist.bin)

### - Jaccard dissimilarity matrix of metabolites 
# only compares shared metabolite presences, ignores shared absences
sampsByCompounds[1:5,1:5]
sampsByCompounds_bin <- (sampsByCompounds > 0) * 1
sampsByCompounds_bin[1:5,1:5]

jaccard<- as.data.frame(as.matrix(vegdist(sampsByCompounds_bin, 
                                          method = "jaccard",diag = T,upper = T,binary = T)))


### ---save once more with the bray curtis distance matrices
save(sampsByCompounds,  sampsByCompounds_bin, sampsCompsStand, diags, cscs, heat.class, heat.itol, 
     braydist,braydist.bin,jaccard,
     file = outfile)



#### END HERE ##########################################################################################
########################################################################################################
########################################################################################################
########################################################################################################