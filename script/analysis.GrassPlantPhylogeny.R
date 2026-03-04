# Analyses to answer Q1: How does the endophyte community vary by host plant?
# Written by Liz Bowman, PhD eabowman@utexas.edu
# Sep. 20, 2024

# Host data ----
host.tree <- read.tree("data/G045_Concatenated_ITStrnL_withbranchlengths.nwk")
data.host.sp <- read.csv('data/G045.Field.data.csv', header = T)

## Calculate host phylogenetic distance ----
# Round the node labels to the tenths
host.tree$node.label <- as.character(round(as.numeric(host.tree$node.label),
                                           digits = 2))

# Get group data for the tips
select(data.host.sp, SampleID, Native.status, Species_final_nov24) -> tree.groups

# Change tip labels
data.frame(SampleID = host.tree$tip.label,
           species = NA) -> new.tip

for(i in 1:nrow(new.tip)){
  sample.id.i <- new.tip[i, 'SampleID']
  as.character(tree.groups$SampleID) -> tree.groups$SampleID
  species.i <- tree.groups[tree.groups$SampleID == sample.id.i, 'Species_final_nov24']
  new.tip[i, 'species'] <- species.i
}

host.tree$tip.label <- new.tip$species

# Plot
ggtree(host.tree) +
  geom_tiplab(size = 3) +
  geom_nodelab(hjust = 1.4, vjust = -0.6, size = 3) +
  theme_tree() +
  # labs(title = "Host Plant: trnL and ITS") +
  xlim(-10, 20)

ggsave('figures/SupplementaryFigureS1.jpeg', plot = last_plot(),
       width = 6, height = 8, units = 'in', device = 'jpeg')
