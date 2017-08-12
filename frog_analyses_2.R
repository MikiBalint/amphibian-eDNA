library(vegan)
library(tidyverse)
library(boral)
library(mvabund)
library(unmarked)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggsn)
library(rgdal)

rm(list=ls())

# Read abundance tables
OwnAssign = read.csv(file="abundances/frogs_16S_own.tab", sep="\t",
                     header=T, row.names=1)
EmblAssign = read.csv(file="abundances/frogs_16S_EMBL.tab", sep="\t",
                      header=T, row.names=1)

# Keep only 16S sequences
EmblAssign = EmblAssign[EmblAssign[,"experiment"] == "16S",]

# Sequence variant metadata
MetaOwn = data.frame(best_id = OwnAssign$best_identity.16S_fasta_nonredundant_noprimers,
                best_match = OwnAssign$best_match, 
                count = OwnAssign$count, 
                family = OwnAssign$family_name,
                genus = OwnAssign$genus_name, 
                sci_name = OwnAssign$scientific_name, 
                taxid = OwnAssign$taxid, 
                sequence = OwnAssign$sequence,
                seq_names = rownames(OwnAssign), 
                row.names=9)
MetaEmbl = data.frame(best_id = EmblAssign$best_identity.db_16S,
                      best_match = EmblAssign$best_match, 
                      count = EmblAssign$count, 
                      family = EmblAssign$family_name,
                      genus = EmblAssign$genus_name, 
                      sci_name = EmblAssign$scientific_name, 
                      taxid = EmblAssign$taxid, 
                      sequence = EmblAssign$sequence,
                      seq_names = rownames(EmblAssign), 
                      row.names=9)
MetaAll = rbind(MetaOwn, MetaEmbl)

# there were taxonomic changes since the bioinformatics
# correct the names
# if the script is adapted to another system: I strongly 
# recommend to cross-check and finalize the names with a taxonomist 
# at this step before moving on: levels(MetaAll$sci_name)
correct_names <- MetaAll$sci_name

correct_names <- gsub("Pseudopaludicola mystacalis", "Pseudopaludicola sp.", 
                      correct_names)
correct_names <- gsub("Dendropsophus leucophyllatus", "Dendropsophus arndti", 
                      correct_names)
correct_names <- gsub("Pseudis limellum", "Lysapsus limellum", 
                      correct_names)
correct_names <- gsub("Phyllomedusa azurea", "Pithecopus azureus", 
                      correct_names)
correct_names <- gsub("Scinax sp. FB-2014a", "Scinax madeirae", 
                      correct_names)
correct_names <- gsub("Dendropsophus minutus", "Dendropsophus minutus A", 
                      correct_names)
correct_names <- gsub("Dendropsophus nanus", "Dendropsophus nanus A", 
                      correct_names)
correct_names <- gsub("Elachistocleis sp.", "Elachistocleis sp. A", 
                      correct_names)
correct_names <- gsub("Hypsiboas punctatus", "Hypsiboas punctatus A", 
                      correct_names)
correct_names <- gsub("Leptodactylus elenae", "Leptodactylus elenae A", 
                      correct_names)
correct_names <- gsub("Physalaemus albonotatus", "Physalaemus cf. albonotatus", 
                      correct_names)
correct_names <- gsub("Physalaemus centralis", "Physalaemus centralis A", 
                      correct_names)
correct_names <- gsub("Scinax ruber", "Scinax cf. nasicus", 
                      correct_names)
correct_names <- gsub("Scinax nasicus", "Scinax nasicus A", 
                      correct_names)
correct_names <- gsub("Pseudis laevis", "Lysapsus laevis", 
                      correct_names)

# taxonomically correct scientific names in the metadata
MetaAll$sci_name <- as.factor(correct_names)

# Obiclean status files
StatusOwn = OwnAssign[, grepl("obiclean.status", colnames(OwnAssign))]
StatusEmbl = EmblAssign[, grepl("obiclean.status", names(EmblAssign))]
  
# Keep 16S statuses
StatusEmbl = (StatusEmbl[, grepl("16", colnames(StatusEmbl))])
StatusAll = rbind(StatusOwn, StatusEmbl)

# status counts
MetaAll = data.frame(MetaAll, 
                     h_count = rowSums(StatusAll == "h", na.rm=T),
                     s_count = rowSums(StatusAll == "s", na.rm=T), 
                     i_count = rowSums(StatusAll == "i", na.rm=T))

# Abundances matrices
AbundOwn = OwnAssign[, grepl("sample", colnames(OwnAssign))]
AbundEmbl = EmblAssign[, grepl("sample", names(EmblAssign))]

# Keep 16S abundances
AbundEmbl = (AbundEmbl[, grepl("16", colnames(AbundEmbl))])
AbundAll = rbind(AbundOwn, AbundEmbl)

# Sequnce variant abundances at least once head
AbundHead = AbundAll[MetaAll$h_count > 0,]
MetaHead = MetaAll[MetaAll$h_count > 0,]

# site metadata
site_meta <- read.csv(file="site_names.csv", row.names = 1)

# NMDS to compare sample replicates and controls
replicate_nmds = metaMDS(t(AbundHead), k=3)

replicate_coordinates <- data.frame(replicate_nmds$points)

# connecting centroids from here:
# https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot
# create a dataframe for calculating centroids
for_centroids <- data.frame(replicate_coordinates[,c(1,2)], 
                            type = site_meta$type)

# re-order df alphabetically according to the type
for_centroids <- for_centroids[order(for_centroids$type),]

# calculate centroid coordinates
replicate_centroids <- aggregate(cbind(MDS1, MDS2) ~ type,
                                 for_centroids,
                                 mean)

replicate_coordinates_centroids <- merge(for_centroids, 
                                         replicate_centroids,
                                         by = "type")
rownames(replicate_coordinates_centroids) <- rownames(for_centroids)

# plot replicates and connect centroids
ggplot(data = replicate_coordinates_centroids) +
  geom_point(mapping = aes(x=MDS1.x, 
                           y=MDS2.x, 
                           color = type)) +
  geom_point(data = replicate_centroids, 
             mapping = aes(x=MDS1,
                           y=MDS2,
                           color = type),
             size=3,
             shape = 1) +
  geom_segment(aes(x=MDS1.y, 
                   y=MDS2.y, 
                   xend=MDS1.x, 
                   yend=MDS2.x,
                   color = type)) +
  geom_text(mapping = aes(x=MDS1.x, 
                          y=MDS2.x, 
                          color = type,
                          label = type),
            nudge_y = 0.14)

# select "weird" replicates: 
# they are at least n times further away from their pond centroid
# than the average distance of all replicates
replicate_centroid_distances <- c()
for (i in rownames(replicate_coordinates_centroids)) {
  single_distance <- sqrt((replicate_coordinates_centroids[i,"MDS1.x"] - 
                             replicate_coordinates_centroids[i,"MDS1.y"])^2 +
                            (replicate_coordinates_centroids[i,"MDS2.x"] - 
                               replicate_coordinates_centroids[i,"MDS2.y"])^2)
  replicate_centroid_distances <- rbind(replicate_centroid_distances, 
                                        single_distance)
}
replicate_coordinates_centroids <- cbind(replicate_coordinates_centroids,
                                         distance = replicate_centroid_distances)

# average distances for sample types
distance_centroids <- aggregate(distance ~ type,
                                replicate_coordinates_centroids,
                                mean)

# include distances into the replicate coordinate object
replicate_coordinates_centroids <- merge(distance_centroids, 
                                         replicate_coordinates_centroids,
                                         by = "type")
rownames(replicate_coordinates_centroids) <- rownames(for_centroids)

# filter replicates by distance from centroid criterium
distance_criterium = 2

distance_filter <- replicate_coordinates_centroids$distance.y < 
  distance_criterium * replicate_coordinates_centroids$distance.x

# ordination plot without the weird replicates
ggplot(data = replicate_coordinates_centroids[distance_filter,]) +
  geom_point(mapping = aes(x=MDS1.x, 
                           y=MDS2.x, 
                           color = type)) +
  geom_point(data = replicate_centroids, 
             mapping = aes(x=MDS1,
                           y=MDS2,
                           color = type),
             size=3,
             shape = 1) +
  geom_segment(aes(x=MDS1.y, 
                   y=MDS2.y, 
                   xend=MDS1.x, 
                   yend=MDS2.x,
                   color = type))

###
# clean the negative controls: remove the maximum number of reads obseved from any negative control
# from the read counts of a sequence variant in each replicate

# get the negative control indices according to the read abundance matrix

# field tap water negative control
CEN <- grep("sample.Cen", names(AbundHead))
colnames(AbundHead)[CEN]

# glass fiber filter extraction control
GEC = grep("sample.NTC", names(AbundHead))
colnames(AbundHead)[GEC]

# nylon filter extraction control
NEC = grep("sample.NC", names(AbundHead))
colnames(AbundHead)[NEC]

# PCR negative controls
NPC = grep("sample.P.NC", names(AbundHead))
colnames(AbundHead)[NPC]

# multiplexing control
MPX = grep("sample.MPX", names(AbundHead))
colnames(AbundHead)[MPX]

# Maximum number of reads in any control sample
MaxControl = apply(AbundHead[,c(NPC,NEC,GEC,MPX, CEN)], 1, max)

# Extract the highest read numbers in a control from every sample
# checks
# AbundHead[1:5,100:102]
# MaxControl[1:5]
# sweep(AbundHead[1:5,100:102], 1, MaxControl[1:5], "-")
AbundControlled = sweep(AbundHead, 1, MaxControl, "-")

AbundControlled[AbundControlled < 0] <- 0

summary(apply(AbundHead,1,sum))

# keep only non-weird replicates (based on the ordination above)
AbundNoWeird <- AbundControlled[,distance_filter]

# keep only replicates with reads
AbundNoWeird <- AbundNoWeird[,apply(AbundNoWeird,2,sum) > 0]

# combine the replicates of samples
# get sample names, code from here: http://stackoverflow.com/questions/9704213/r-remove-part-of-string
SampleNames = levels(as.factor(sapply(strsplit(names(AbundNoWeird), 
                                               split='16S', fixed=TRUE), 
                                      function(x) (x[1]))))

# Sum the replicates for each sample
SummedReps = data.frame(row.names = rownames(AbundNoWeird))
for (i in 1:length(SampleNames)){
  ActualSet = grep(SampleNames[i], names(AbundNoWeird)) # grep the columns of interest
  SummedReps = cbind(SummedReps, apply(AbundNoWeird[ActualSet], 1, sum))
}
colnames(SummedReps) = SampleNames

###
# false positives: the presence of a sequence variant is accepted if 
# found in at least two replicates

# In how many replicates observed per sample?
PresentReps = data.frame(row.names = rownames(AbundNoWeird))
for (i in 1:length(SampleNames)){
  ActualSet = grep(SampleNames[i], names(AbundNoWeird))
  Selected = AbundNoWeird[ActualSet]
  Selected[Selected > 0] <- 1 # set the read numbers to 1
  PresentReps = cbind(PresentReps, apply(Selected, 1, sum))
}
colnames(PresentReps) = SampleNames

# Set read numbers to 0 in a sample if the sequence variant was not observed in at least
# two PCR replicates
SummedControlled = SummedReps
SummedControlled[PresentReps < 2] <- 0

# remove empty sequence variants
SummedControlled <- SummedControlled[apply(SummedControlled,1,sum) > 0,]

# adjust the sequence variant metadata
MetaControlled <- MetaHead[rownames(MetaHead) %in% rownames(SummedControlled),]

# Categories for sequence variants: frogs, humans, other aquatic stuff, etc.
# Amphi, HighGroup, "Homo sapiens", FarmAnim, Bird, Fish, Insect, Mammal
# from the MetaHead$sci_name variable

Amphi = c("Dendropsophus arndti","Dendropsophus melanargyreus",
          "Dendropsophus minutus A","Dendropsophus nanus A","Dendropsophus salli","Dermatonotus muelleri",
          "Elachistocleis sp. A","Eupemphix nattereri","Hylinae","Hyloidea",
          "Hypsiboas geographicus",
          "Hypsiboas punctatus A","Hypsiboas raniceps","Leptodactylus",
          "Leptodactylus chaquensis","Leptodactylus elenae A","Leptodactylus fuscus",
          "Leptodactylus podicipinus","Leptodactylus syphax","Leptodactylus vastus",
          "Osteocephalus","Osteocephalus taurinus","Pithecopus azureus",
          "Phyllomedusa boliviana","Physalaemus cf. albonotatus","Physalaemus centralis A",
          "Physalaemus cuvieri","Lysapsus limellum","Pseudis paradoxa",
          "Pseudopaludicola sp.","Rhinella schneideri","Scinax fuscomarginatus",
          "Scinax fuscovarius","Scinax nasicus A","Scinax cf. nasicus","Scinax madeirae",
          "Bufonidae","Elachistocleis","Gastrophryninae","Hylidae",
          "Leptodactylus latinasus","Melanophryniscus","Pseudis","Lysapsus laevis","Scinax")

HighGroup = c("root", "Bilateria", "Amniota")

FarmAnim = c("Canis", "Sus scrofa", "Gallus", "Mus", "Bos", "Capreolus capreolus")

Bird = c("Meleagris gallopavo", "Jacana", "Gallus")

Fish = c("Astyanax","Brachyplatystoma rousseauxii","Callichthys callichthys","Characidae",
         "Characiphysae","Opisthonema oglinum","Salmonidae","Sardinops","Sardinops sagax",
         "Siluriformes","Siluroidei","Synbranchus marmoratus","Thraupidae","Thunnini")

Insect = c("Micrathyria ocellata", "Libellulidae", "Micrathyria", "Tabanus", "Tramea")

Mammal = c("Homo sapiens", "Canis", "Sus scrofa", "Mus", "Bos", "Capreolus capreolus",
           "Myotis", "Laurasiatheria")

# Proportions of different groups
pdf("Fig_pie.pdf", width = 3.15, height = 3.15)
par(mar=c(1,1,1,1))
pie(c(sum(SummedControlled[MetaControlled$sci_name %in% Amphi,]), 
      sum(SummedControlled[MetaControlled$sci_name %in% Mammal,]),
      sum(SummedControlled[MetaControlled$sci_name %in% Fish,]),
      sum(SummedControlled[MetaControlled$sci_name %in% Insect,]),
      sum(SummedControlled[MetaControlled$sci_name %in% Bird,]),
      sum(SummedControlled[MetaControlled$sci_name %in% HighGroup,])),
    labels = paste(c("Frogs\n2 158 534", "Mammals\n14 006", 
                     "Fish\n1 692 613", "Insects\n304 059", 
                     "Birds\n967", "Higher groups\n1 063 156")),
    col = c(gray(0.9), gray(0.75), gray(0.60), gray(0.45), gray(0.30), gray(0.15)),
    cex = 0.8)
dev.off()

# Farm animals:
sum(SummedControlled[MetaHead$sci_name %in% FarmAnim,])

# Humans:
sum(SummedControlled[MetaHead$sci_name == "Homo sapiens",])

# Libellulidae: 
sum(SummedControlled[MetaHead$sci_name %in% c("Micrathyria ocellata", "Libellulidae", 
                                            "Micrathyria", "Tramea"),])

# Aggregate frog sequence variants according to the species 
FrogAggregate = data.frame(name = MetaControlled$sci_name[MetaControlled$sci_name %in% Amphi], 
                         SummedControlled[MetaControlled$sci_name %in% Amphi,])
FrogCounts = aggregate(. ~ name, FrogAggregate, sum, na.action = na.exclude)
rownames(FrogCounts) = FrogCounts$name
FrogCounts = FrogCounts[,2:79]

# Remove ambiguous assignments:
ambiguous = c("Hylinae", "Hylidae", "Hyloidea", "Leptodactylus", "Osteocephalus",
              "Bufonidae", "Gastrophryninae", "Pseudis", "Scinax")

# Species per sample dataset without ambiguous identifications
FrogCounts = FrogCounts[!(rownames(FrogCounts) %in% ambiguous),]

# Remove taxa with no observations 
# (can happen due to the cleanup steps on the sequence variants)
FrogCounts = FrogCounts[apply(FrogCounts, 1, sum) != 0, ]

# Generate site metadata from column names
# Lake codes
ponds = substring(names(FrogCounts), 8, 9)

# Pond metadata: ponds, preservation/extraction methods, read numbers
PondMeta = data.frame(ponds = ponds, 
                      type = c(rep(c("control"),15), 
                               rep(c("A","B","C"),12),
                                  rep("A",4), rep("B",3), rep("C",3),
                                  rep("A",3), rep("B",3), rep("C",3),
                                  rep("A",3), rep("B",3), rep("C",2)),
                      reads = apply(FrogCounts, 2, sum))
PondCodes = c("T1", "T2", "T3", "T4", "T5")

###
# Establish a sequence rarity threshold with the positive controls

# Get positive samples
PositiveCodes = c("PCE", "PCUNE", "Lvastus")

# Grep names of positive and field tap water samples. 
# Sample names are in the 'lakeCodes'
IsPosField = grep(paste(PositiveCodes, collapse='|'), 
                  names(FrogCounts), ignore.case=TRUE)

# Create read abundance table for the positive controls
FrogsInPosField = FrogCounts[,IsPosField]
FrogsInPosField = FrogsInPosField[apply(FrogsInPosField,1,sum) > 0,]

# list of species included into the positives
# Species in the positiv control
PosList = sort(c("Pithecopus azureus",
                 "Physalaemus centralis A",
                 "Hypsiboas geographicus",
                 "Pseudopaludicola sp.",
                 "Hypsiboas punctatus A",
                 "Leptodactylus fuscus",
                 "Hypsiboas raniceps",
                 "Lysapsus limellum",
                 "Leptodactylus podicipinus",
                 "Scinax madeirae",
                 "Dendropsophus arndti",
                 "Dendropsophus nanus A",
                 "Leptodactylus vastus")) # L. vastus was a field positive control: present in an aquarium from which water was used as field positive.

# comparison list with read numbers to find a false positive threshold
FrogsInPosField <- data.frame(FrogsInPosField,
                              max_positives = apply(FrogsInPosField, 1, max),
                              max_positives = apply(FrogsInPosField, 1, max),
                              is_positive = rownames(FrogsInPosField) %in% PosList)

# the maximum read number of a non-positive control species in a positive is 32 
# the false positive threshold is then 40.
# Physalaemus cuvieri should not be in the positives. It has read numbers 
# that are typical for true positives (e.g. Hypsiboas raniceps), so it is not 
# used as a positive threshold.
positive_threshold <- 40

FrogCounts[FrogCounts < 40] <- 0
FrogCounts <- FrogCounts[apply(FrogCounts,1,sum) > 0,]
FrogCounts <- FrogCounts[,apply(FrogCounts,2,sum) > 0]

###
# summed frog read abundances in ponds

# Grep names of all samples. Sample names are in the 'PondCodes'
IsPondSample = grep(paste(PondCodes, collapse='|'), 
                    names(FrogCounts), ignore.case=TRUE)

# Create read abundance table for frog species found in the ponds
FrogsInPonds = FrogCounts[,IsPondSample]
FrogsInPonds = FrogsInPonds[apply(FrogsInPonds,1,sum) > 0,]

# Sum up species counts by lakes
FrogPonds = data.frame(row.names = rownames(FrogsInPonds))
for (i in 1:length(PondCodes)) {
  ActualSet = grep(PondCodes[i], names(FrogsInPonds))
  FrogPonds = cbind(FrogPonds, apply(FrogsInPonds[ActualSet], 1, sum))
}
colnames(FrogPonds) = PondCodes

# summed read abundances of species in the five ponds
write.csv(file="Frog_sums_ponds.csv", 
          cbind(FrogPonds, Total = apply(FrogPonds,1,sum)))

#####
# Ecological signal: differences among the lakes VS preservation/extraction methods
# Transpose the frog counts: species in columns, samples in lines
frogs_in_ponds <- t(FrogsInPonds)

# Filter for lake samples
pond_filter <- rownames(PondMeta) %in% rownames(frogs_in_ponds)
pond_only_meta <- PondMeta[pond_filter,]

# Model-based comparison of the lakes
frogs_mvabund = mvabund(frogs_in_ponds)

# Multispecies GLM and ANOVA. Reads control for differences in sequencing depth.
manyglm_m1 = manyglm(frogs_mvabund ~ ponds, 
             data=pond_only_meta,
             family = "negative.binomial")

# analysis of variance, pond effect
m1.anova = anova(manyglm_m1, nBoot = 100, p.uni = "adjusted")

# nice anova table
m1.anova$table

# Visualize the community results
# Model-based ordination with reads as covariates

# Overdispersion parameter
summary(manyglm_m1$theta)
hist(manyglm_m1$theta)

# overdispersion prior
set.prior <- list(type = c("normal","normal","normal","uniform"),
                  hypparams = c(100, 20, 100, 3))

# model-based ordination
boral_m1 <- boral(frogs_in_ponds,
                  family = "negative.binomial", 
                  num.lv = 2,
                  prior.control = set.prior)

# group centroids of lakes
LV_lake_factors <- data.frame(boral_m1$lv.median, 
                              ponds=factor(pond_only_meta$ponds))
pond_centroids <- data.frame()
for (i in c("T1", "T2", "T3", "T4", "T5")){
  LV1_centroid <- mean(LV_lake_factors$lv1[LV_lake_factors$ponds == i])
  LV2_centroid <- mean(LV_lake_factors$lv2[LV_lake_factors$ponds == i])
  pond_centroids <- rbind(pond_centroids, c(LV1_centroid, LV2_centroid))
}
colnames(pond_centroids) = c("LV1_centroid", "LV2_centroid")
rownames(pond_centroids) <- c("T1", "T2", "T3", "T4", "T5")

### VAES ordination
# pond ordination with VAES presence-absences
VAES = read.csv(file = "abundances/VAES_presence-absence.csv", 
                header = T, row.names = 1)

# # VAES with LVM ordiantion
boral_VAES_m1 <- boral(t(VAES),
                       family = "binomial",
                       num.lv = 2)

# Colors for ordination
palette(colors())
LakeCol = gsub("T1", "green", pond_only_meta$ponds)
LakeCol = gsub("T2", "orange", LakeCol)
LakeCol = gsub("T3", "red", LakeCol)
LakeCol = gsub("T4", "purple", LakeCol)
LakeCol = gsub("T5", "blue", LakeCol)
MyColors = data.frame(cols = c("green","orange","red","purple","blue"),
                      ponds = levels(factor(pond_only_meta$ponds)))

# ordination plots
pdf("ordinations.pdf", width = 3.15, height = 6)
par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(1,1,0,0))
# layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE),
#        heights=c(2,1))
ordifrog <- ordiplot(boral_m1$lv.median, type = "none", cex =2, 
                     display = "sites", # xlim = c(-0.1,0.1), ylim = c(-6,7),
                     cex.lab = 1.2, cex.axis = 1.2,
                     xlab = "Community axis 1",
                     ylab = "Community axis 2")
points(ordifrog,"sites", pch=20 ,col=LakeCol, cex = 1.5)
ordiellipse(ordifrog, factor(pond_only_meta$ponds), cex=.5, 
            draw="polygon", col="green",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T1")))
ordiellipse(ordifrog, factor(pond_only_meta$ponds), cex=.5, 
            draw="polygon", col="orange",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T2")))
ordiellipse(ordifrog, factor(pond_only_meta$ponds), cex=.5, 
            draw="polygon", col="red",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T3")))
ordiellipse(ordifrog, factor(pond_only_meta$ponds), cex=.5, 
            draw="polygon", col="purple",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T4")))
ordiellipse(ordifrog, factor(pond_only_meta$ponds), cex=.5, 
            draw="polygon", col="blue",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T5")))
points(pond_centroids, pch = 3)
text(pond_centroids, labels = rownames(pond_centroids),
     pos = 4, offset = 0.3, cex = 0.8)
plot(procrustes(pond_centroids, boral_VAES_m1$lv.median), 
     main = "", kind = 0,
     xlab = "Community axis 1",
     ylab = "Community axis 2",
     cex.lab = 1.2, cex.axis = 1.2)
points(procrustes(pond_centroids, boral_VAES_m1$lv.median), 
       pch = 21, bg = as.character(MyColors$cols), cex = 1.5)
lines(procrustes(pond_centroids, boral_VAES_m1$lv.median), 
       type = c("arrows"), angle = 20, length = 0.1)
text(procrustes(pond_centroids, boral_VAES_m1$lv.median),
     pos = 4, offset = 0.5, cex = 0.8)
dev.off()

# permuation test of the correspondence
protest(pond_centroids, boral_VAES_m1$lv.median)

###
# recovery of species in the positive controls
# PCE: DNA from 12 species in equal concentrations
C.PCE = grep("sample.P.PCE", names(FrogCounts))

# PCUNE: DNA from 12 species in stepwise doubled concentrations
C.PCUNE = grep("sample.P.PCUNE", names(FrogCounts))

# Species in the positiv control
PosList_mock = sort(c("Pithecopus azureus",
                      "Physalaemus centralis",
                      "Hypsiboas geographicus",
                 "Pseudopaludicola sp.",
                 "Hypsiboas punctatus",
                 "Leptodactylus fuscus",
                 "Hypsiboas raniceps",
                 "Lysapsus limellum",
                 "Leptodactylus podicipinus",
                 "Scinax madeirae",
                 "Dendropsophus arndti",
                 "Dendropsophus nanus"))

# non-equimolar positive control DNA concentrations
PCEUNEConc = read.csv(file="PCUNE_conc.csv", header=T, row.names = 1)
PCEUNEConc <- PCEUNEConc[order(rownames(PCEUNEConc)),]
PCEUNEConc <- data.frame(PCEUNEConc, FrogCounts[PosList_mock, C.PCUNE])
PCEUNEConc <- PCEUNEConc[order(PCEUNEConc$conc),]

# The seven PCE samples are positively correlated at R > 0.7.
hist(cor(t(FrogCounts[PosList_mock,C.PCE])), 
     xlab = "Correlation coefficient",
     main = "")

# plot positives
# ranges for non-equimolar positive plot
RangeX = range(PCEUNEConc[,"conc"])
RangeY = range(apply(FrogCounts[,C.PCUNE],1,max))

# reorder equimolar positive
PosList_mock[order(rownames(PCEUNEConc))]

# pdf("positives.pdf", width = 5, height = 8)
palette(colors())
par(mfrow = c(2,1), mar=c(6,5,1,1))
plot(c(1:12), seq(0.1, (max(FrogCounts[,C.PCE])), 
                  (max(FrogCounts[,C.PCE]))/12), type="n",
     xaxt="n", 
     xlab="", ylab="Read numbers") 
axis(1, at = c(1:12), labels = rownames(PCEUNEConc), las = 2, cex.axis=0.6)
for (i in C.PCE) {
  lines(c(1:12), (FrogCounts[rownames(PCEUNEConc),i]), 
        col=i*10, lwd=2, type = "o", pch = 19)
}
par(mar=c(5,5,3,1))
plot(RangeX, RangeY, type="n", 
     xlab="DNA concentrations (ng/ul)", ylab = "Read numbers", 
     log = "x")
for (i in 3:9){
  lines(na.omit(data.frame(PCEUNEConc$conc, PCEUNEConc[,i])),
        lwd=2, col = i*10, type = "o", pch = 19)
}

###
# Species occupancy model
#### SCRIPT TO RUN THE FALSE POSITIVE SPECIES OCCUPANCY MODEL ####
#### USED TO ESTIMATE THE DETECTION AND FALSE POSITIVE PROBABILITY OF EACH EDNA METHOD ####
#### SOM script by Thierry Chambert
options(scipen=999)

### 1. EXTRACT AND FORMAT THE DATA
# load csv files
dat0 <- read.csv("som_data_input.csv")
dat0$Site = as.numeric( gsub("T", "", dat0$Site) )
dat0$eDNA.Method = as.character(dat0$eDNA.Method)

## Spatial Replicates
for(a in unique(dat0$Site)){
  repl = dat0$Spatial.Replicate[dat0$Site == a]
  mn = min( dat0$Spatial.Replicate[dat0$Site == a] ) - 1
  dat0$Spatial.Replicate[dat0$Site == a] <- repl-mn
}

for(a in unique(dat0$Site)){
  repl = dat0$Spatial.Replicate[dat0$Site == a]
  mn = min( dat0$Spatial.Replicate[dat0$Site == a] ) - 1
  dat0$Spatial.Replicate[dat0$Site == a] <- repl-mn
} # a

for(a in 1:nrow(dat0)){
  splt = strsplit(x=dat0$eDNA.Method[a], split="")
  dat0$eDNA.Method[a] <- splt[[1]][1]
  if( length(splt[[1]]) > 1)
    # dat0$Spatial.Replicate[a] <- ( dat0$Spatial.Replicate[a] * as.numeric(splt[[1]][2]) ) + ( (as.numeric(splt[[1]][2])-1)*10 )
    dat0$Spatial.Replicate[a] <- ( dat0$Spatial.Replicate[a] * as.numeric(splt[[1]][2]) ) + ( (dat0$Spatial.Replicate[a]-1)*10 )
} # a

rk = cbind( rank( unique(dat0$Spatial.Replicate) ), unique(dat0$Spatial.Replicate) )
rk = rk[order(rk[,1]),]
for(a in 1:nrow(dat0)) dat0$Spatial.Replicate[a] <- rk[(rk[,2] == dat0$Spatial.Replicate[a]),1]

## ORDER DATA
dat1 = dat0[ order(dat0$Species, dat0$Site, dat0$eDNA.Method, dat0$Spatial.Replicate, dat0$PCR.Replicate), ]
dat2 = dat1[,c(2:7)]

## eDNA Method as numeric
dat2$eDNA.Method[dat2$eDNA.Method == "A"] <- 1
dat2$eDNA.Method[dat2$eDNA.Method == "B"] <- 2
dat2$eDNA.Method[dat2$eDNA.Method == "C"] <- 3
dat2$eDNA.Method <- as.numeric(dat2$eDNA.Method)

## Dimension DATA MATRIX
# Columns
UU = unique( cbind(dat2$eDNA.Method, dat2$Spatial.Replicate, dat2$PCR.Replicate) )
U <- UU[order(UU[,1],UU[,2],UU[,3]),]
NCOL = nrow(U)

# Rows
UROW = unique( cbind(dat2$Species, dat2$Site) )
NROW = nrow(UROW)
length(unique(dat2$Species)) * length(unique(dat2$Site))

## Create DATA MATRIX
som.dat = matrix(NA,NROW,NCOL)
colnames(som.dat) = paste(U[,1],U[,2],U[,3],sep=".")

## FILL IN DATA
for(j in 1:NROW){
  spec = levels(dat2$Species)[ UROW[j,1] ]
  site = UROW[j,2]
  c1 = dat2[dat2$Species == spec & dat2$Site == site, "eDNA.Method" ]
  c2 = dat2[dat2$Species == spec & dat2$Site == site, "Spatial.Replicate" ]
  c3 = dat2[dat2$Species == spec & dat2$Site == site, "PCR.Replicate" ]
  cname = paste(c1,c2,c3,sep=".")
  cdat = dat2[as.character(dat2$Species) == spec & dat2$Site == site, "Detection" ]
  fill = which(colnames(som.dat) %in% cname)
  som.dat[j,fill] <- cdat
}# j

### eDNA methods ###
eDNA.Meth = list()
# Method A
eDNA.Meth[[1]] <- which(U[,1] == 1)
# Method B
eDNA.Meth[[2]] <- which(U[,1] == 2)
# Method C
eDNA.Meth[[3]] <- which(U[,1] == 3)

### OTHER COVARIATES: Species ID
sp.cod = as.numeric(dat2$Species)
sp.cod[sp.cod < 10 ] <- paste("0", sp.cod[sp.cod < 10 ], sep="")
info = data.frame(unique( cbind(as.character(dat2$Species), as.character(sp.cod), as.numeric(dat2$Site)) ), stringsAsFactors = F )
names(info) = c("species.name","species.code", "site") 

### Distinguish FP/TP detections
## IDENTIFY DETECTION DATA by PCR 
pcr.id = NA
nb.pcr = ncol(som.dat) / 4
for(jj in 1:nb.pcr) pcr.id = c( pcr.id , rep(jj, 4) )
pcr.id = pcr.id[-1]

## Sum detection over each PCR
pcr.sum = matrix(NA,nrow(som.dat),nb.pcr)
for(a in 1:nb.pcr) pcr.sum[,a] = apply(som.dat[,pcr.id==a], 1, sum, na.rm=T)

## LUMP the 4 PCR replicate as ONE DETECTION DATA
## AND CODE UNAMBIGUOUS DETECTION AS "2" // "1" is for AMBIGUOUS DETECTIONS
data.1 = pcr.sum
data.1[data.1 >= 2] <- 2


#####################
#### 2. ANALYSIS ####
#####################

## COVARIATES 
# Detection BY eDNA METHOD
eDNA.Meth.all = eDNA.Meth
eDNA.Meth[[1]] <- seq(1, max(eDNA.Meth.all[[1]])/4 ,1)
eDNA.Meth[[2]] <- seq(max(eDNA.Meth[[1]])+1, max(eDNA.Meth.all[[2]])/4 ,1)
eDNA.Meth[[3]] <- seq(max(eDNA.Meth[[2]])+1, max(eDNA.Meth.all[[3]])/4 ,1)

METH1 <- METH2 <- METH3 <- matrix(NA, nrow(data.1),ncol(data.1))
METH1[,eDNA.Meth[[1]]] = 1
METH1[,eDNA.Meth[[2]]] = 0
METH1[,eDNA.Meth[[3]]] = 0
METH2[,eDNA.Meth[[1]]] = 0
METH2[,eDNA.Meth[[2]]] = 1
METH2[,eDNA.Meth[[3]]] = 0
METH3[,eDNA.Meth[[1]]] = 0
METH3[,eDNA.Meth[[2]]] = 0
METH3[,eDNA.Meth[[3]]] = 1
occ <- list(meth1 = METH1, meth2 = METH2, meth3 = METH3)

# Separated "sites" BY SPECIES 
site.cov <- data.frame( species = factor(info$species.code) )

## FORMAT DATA for "unmarked"
data <- unmarkedFrameOccuFP(y = data.1, siteCovs=site.cov, obsCovs=occ, type = c(0,0,3))


## BUILD/RUN MODEL
fm1 <- occuFP(detformula = 	~ 0 + meth1 + meth2 + meth3, 
              FPformula = ~ 0 + meth1 + meth2 + meth3, 
              Bformula = 	~ 0 + meth1 + meth2 + meth3,
              stateformula = ~ 0 + species, data = data)


#####################
#### 3. RESUTLS  ####
#####################
## ESTIMATES
# Detection (det)
det = cbind( round(plogis( coef(fm1, type="det") ), 3) ,
             round( plogis( confint(fm1, type="det") ), 3) )
colnames(det) = c("estimate", "2.5%CI", "97.5%CI")
det

# FP rate
fp = cbind( round(plogis( coef(fm1, type="fp") ), 3) ,
            round( plogis( confint(fm1, type="fp") ), 3) )
colnames(fp) = c("estimate", "2.5%CI", "97.5%CI")
fp

## plot the SOM results
# Coefficient plots
par(mfrow=c(2,1), mar=c(3,10,2,0.5))
plot(c(min(det), mean(det), max(det)), c(1:3),
     yaxt = "n",
     ylab = "", xlab = "", main = "Detection probability",
     type = "n",
     ylim = c(0.5,3.5))
points(det[,"estimate"], c(nrow(det):1),
       pch = 19, cex = 1.2)
segments(det[,"2.5%CI"], c(nrow(det):1), 
         det[,"97.5%CI"], c(nrow(det):1),
         lwd = 1.2)
axis(2, at=c(nrow(det):1), 
     label=c("GFF (2 μm, in liquid, A)",
             "GFF (2 μm, dried: B)",
             "GFF (0.2 μm, dried: C)"),
     las=1, cex.axis=0.9,
     tick=F)

plot(c(min(fp), mean(fp), max(fp)), c(1:3),
     yaxt = "n",
     ylab = "", xlab = "", main = "False positives",
     type = "n",
     ylim = c(0.5,3.5))
points(fp[,"estimate"], c(nrow(fp):1),
       pch = 19, cex = 1.2)
segments(fp[,"2.5%CI"], c(nrow(fp):1), 
         fp[,"97.5%CI"], c(nrow(fp):1),
         lwd = 1.2)
axis(2, at=c(nrow(fp):1), 
     label=c("GFF (2 μm, in liquid, A)",
             "GFF (2 μm, dried: B)",
             "GFF (0.2 μm, dried: C)"),
     las=1, cex.axis=0.9,
     tick=F)

###
# Map
# area box
chiquitos_bbox <- make_bbox(lat = c(-16.38, -16.356), 
                            lon = c(-62.01, -61.995))

# pond coordinates
ponds <- data.frame(row.names = c("T1", "T2", "T3", "T4", "T5"),
                    lon = c(-62.00214, -62.000538, -62.00275, -62.005674, -62.001422),
                    lat = c(-16.3767, -16.357364, -16.3607, -16.366689, -16.359885)
                    )

# google terrain area
chiquitos_goo <- get_googlemap(center = c(-62.002, -16.3665), 
                               source = "google", 
                               maptype = "satellite",
                               size = c(320,640),
                               scale = 2,
                               zoom = 15
                               )

# pond polygons
ponds_polygon <- readOGR("pond_polygons.kml")

# png(filename = "ponds.png", width = 6, height = 8, units = "cm", res = 300)
ggmap(chiquitos_goo) + 
  geom_polygon(data = ponds_polygon, 
               aes(x = long, 
                   y = lat, 
                   group = group),
               color = "white", 
               fill = "lightblue", 
               alpha = 0.7) +
  labs(x = "Longitude",
     y = "Latitude") +
  geom_text(data = ponds, 
            aes(label = paste("  ", as.character(rownames(ponds)), sep="")), 
            hjust = 0.2, 
            color = "white",
            size = 3)
# dev.off()
