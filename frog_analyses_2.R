# OBITools commands

# assignment with EMBL only
# ecotag -d /phylodata/mbalint/databases/ecoPCR_embl_130/ecopcr_embl_130 -R -m 0.98 ../09_ecopcr/NCBI/16S/db_16S.fasta frogs_clean.fasta > frogs_16S_only_EMBL_assigned.fasta

library(vegan)
library(tidyverse)
library(boral)
library(mvabund)
library(unmarked)
# # library(vegan3d)
# library(bvenn)
# library(knitr)
# library(boral)

# # library(geosphere)
# library(car)
# library(ape) # installed with ctv, infos here: http://www.phytools.org/eqg/Exercise_3.2/

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
  # geom_text(mapping = aes(x=MDS1.x, 
  #                         y=MDS2.x, 
  #                         color = type,
  #                         label = rownames(replicate_coordinates_centroids)))

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
# write.csv(file="test.csv", AbundControlled)

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

Amphi = c("Dendropsophus arndti","Dendropsophus leucophyllatus","Dendropsophus melanargyreus",
          "Dendropsophus minutus","Dendropsophus nanus","Dendropsophus salli","Dermatonotus muelleri",
          "Elachistocleis sp.","Eupemphix nattereri","Hylinae","Hyloidea",
          "Hypsiboas geographicus",
          "Hypsiboas punctatus","Hypsiboas raniceps","Leptodactylus",
          "Leptodactylus chaquensis","Leptodactylus elenae","Leptodactylus fuscus",
          "Leptodactylus podicipinus","Leptodactylus syphax","Leptodactylus vastus",
          "Osteocephalus","Osteocephalus taurinus","Phyllomedusa azurea",
          "Phyllomedusa boliviana","Physalaemus albonotatus","Physalaemus centralis",
          "Physalaemus cuvieri","Pseudis limellum","Pseudis paradoxa",
          "Pseudopaludicola mystacalis","Rhinella schneideri","Scinax fuscomarginatus",
          "Scinax fuscovarius","Scinax nasicus","Scinax ruber","Scinax sp. FB-2014a",
          "Bufonidae","Elachistocleis","Gastrophryninae","Hylidae",
          "Leptodactylus latinasus","Melanophryniscus","Pseudis","Pseudis laevis","Scinax")

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
# pdf("Fig_pie.pdf")
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
    col = c(gray(0.9), gray(0.75), gray(0.60), gray(0.45), gray(0.30), gray(0.15)))
# dev.off()

sum(SummedControlled[MetaHead$sci_name %in% Amphi,])+ 
sum(SummedControlled[MetaHead$sci_name %in% Mammal,])+
sum(SummedControlled[MetaHead$sci_name %in% Fish,])+
sum(SummedControlled[MetaHead$sci_name %in% Insect,])+
sum(SummedControlled[MetaHead$sci_name %in% Bird,])+
sum(SummedControlled[MetaHead$sci_name %in% HighGroup,])

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
PosList = sort(c("Phyllomedusa azurea",
                 "Physalaemus centralis",
                 "Hypsiboas geographicus",
                 "Pseudopaludicola mystacalis",
                 "Hypsiboas punctatus",
                 "Leptodactylus fuscus",
                 "Hypsiboas raniceps",
                 "Pseudis limellum",
                 "Leptodactylus podicipinus",
                 "Scinax madeirae",
                 "Dendropsophus leucophyllatus",
                 "Dendropsophus nanus",
                 "Leptodactylus vastus"))

# comparison list with read numbers to find a false positive threshold
FrogsInPosField <- data.frame(FrogsInPosField,
                              max_positives = apply(FrogsInPosField, 1, max),
                              sum_positives = apply(FrogsInPosField, 1, sum),
                              is_positive = rownames(FrogsInPosField) %in% PosList)

# Physalaemus cuvieri is present only in positive controls, this is weird! 
# May be a mis-naming of H. geographicus?
sum(FrogCounts["Physalaemus cuvieri",])


# the maximum read number of a non-positive control species in a positive is 32. 
# the false positive threshold is then 40 - a bit higher.
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

# modify the names of the species so it matches the latest nomenclature 
# (Jansen et al. 2011, Duellmann et al. 2016, Caminer et al. 2017)
actual_names <- c("Dendropsophus melanargyreus",
                  "Dendropsophus minutus A",
                  "Dendropsophus nanus A",
                  "Dermatonotus muelleri",
                  "Elachistocleis sp. A",
                  "Eupemphix nattereri",
                  "Hypsiboas punctatus A",
                  "Hypsiboas raniceps",
                  "Leptodactylus elenae A",
                  "Leptodactylus syphax",
                  "Leptodactylus vastus",
                  "Osteocephalus taurinus",
                  "Pithecopus azureus",
                  "Phyllomedusa boliviana",
                  "Physalaemus cf. albonotatus",
                  "Physalaemus centralis A",
                  "Pseudis paradoxa",
                  "Rhinella schneideri",
                  "Scinax fuscomarginatus",
                  "Scinax fuscovarius",
                  "Scinax nasicus A",
                  "Scinax cf. nasicus",
                  "Scinax madeirae",
                  "Leptodactylus latinasus",
                  "Lysapsus laevis"
)

rownames(FrogsInPonds)

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

# filter singleton and doubleton species
single_doubleton_filter <- specnumber(t(frogs_in_ponds)) > 2

# Model-based comparison of the lakes
frogs_mvabund = mvabund(frogs_in_ponds[,single_doubleton_filter])

# Multispecies GLM and ANOVA. Reads control for differences in sequencing depth.
manyglm_m1 = manyglm(frogs_mvabund ~ ponds, 
             data=pond_only_meta,
             family = "negative.binomial")

# analysis of variance, pond effect, increase nBoot to 1000
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

# model-based ordination - only with non-singleton species
boral_m1 <- boral(frogs_in_ponds[,single_doubleton_filter], 
                  family = "negative.binomial", 
                  num.lv = 2,
                  prior.control = set.prior)

# Diagnostics
par(mfrow = c(2,2))
plot(boral_m1)

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

# VAES with LVM ordiantion
boral_VAES_m1 <- boral(t(VAES), 
                       family = "binomial", 
                       num.lv = 2)

# VAES_ord <- metaMDS(t(VAES), distance = "jaccard", trymax = 1000)

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
pdf("ordinations.pdf", width = 5, height = 8)
par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(1,1,0,0))
# layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE),
#        heights=c(2,1))
ordifrog <- ordiplot(boral_m1$lv.median, type = "none", cex =2, 
                     display = "sites", # xlim = c(-0.1,0.1), ylim = c(-6,7),
                     cex.lab = 1.2, cex.axis = 1.2,
                     xlab = "Community axis 1",
                     ylab = "Community axis 2")
points(ordifrog,"sites", pch=20 ,col=LakeCol, cex = 2)
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
     pos = 4, offset = 0.3)
plot(procrustes(pond_centroids, boral_VAES_m1$lv.median), 
     main = "", kind = 0,
     xlab = "Community axis 1",
     ylab = "Community axis 2",
     cex.lab = 1.2, cex.axis = 1.2)
points(procrustes(pond_centroids, boral_VAES_m1$lv.median), 
       pch = 21, bg = as.character(MyColors$cols), cex = 2)
lines(procrustes(pond_centroids, boral_VAES_m1$lv.median), 
       type = c("arrows"), angle = 20, length = 0.1)
text(procrustes(pond_centroids, boral_VAES_m1$lv.median),
     pos = 1, offset = 0.7)
dev.off()

# permuation test of the correspondence
protest(pond_centroids, boral_VAES_m1$lv.median)

# protest(pond_centroids, VAES_ord$points)

###
# recovery of species in the positive controls
# PCE: DNA from 12 species in equal concentrations
C.PCE = grep("sample.P.PCE", names(FrogCounts))

# PCUNE: DNA from 12 species in stepwise doubled concentrations
C.PCUNE = grep("sample.P.PCUNE", names(FrogCounts))

# Species in the positiv control
PosList_mock = sort(c("Phyllomedusa azurea",
                      "Physalaemus centralis",
                      "Hypsiboas geographicus",
                 "Pseudopaludicola mystacalis",
                 "Hypsiboas punctatus",
                 "Leptodactylus fuscus",
                 "Hypsiboas raniceps",
                 "Pseudis limellum",
                 "Leptodactylus podicipinus",
                 "Scinax madeirae",
                 "Dendropsophus leucophyllatus",
                 "Dendropsophus nanus"))

# non-equimolar positive control DNA concentrations
PCEUNEConc = read.csv(file="PCUNE_conc.csv", header=T, row.names = 1)
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

pdf("positives.pdf", width = 5, height = 8)
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
dev.off()

###
# species observation models - SOM

# format data
# filter for pond sample PCR replicates
is_pond_replicate <- grep("sample\\.T", names(AbundControlled))
pond_replicate_abunds <- AbundControlled[,is_pond_replicate]
pond_replicate_abunds <- pond_replicate_abunds[apply(pond_replicate_abunds,1,sum) > 0,]

# metadata for pond replicates
pond_replicate_meta <- MetaHead[rownames(MetaHead) %in% rownames(pond_replicate_abunds),]

# aggregate according to scientific names
aggregate_for_som <- data.frame(name = pond_replicate_meta$sci_name[pond_replicate_meta$sci_name %in% Amphi], 
                                pond_replicate_abunds[pond_replicate_meta$sci_name %in% Amphi,])
species_for_som <- aggregate(. ~ name, aggregate_for_som, sum, na.action = na.exclude)
# rownames(species_for_som) = species_for_som$name
# species_for_som = species_for_som[,2:ncol(species_for_som)]

# keep only reliable species - those that were kept after all cleanup steps
species_for_som <- species_for_som[species_for_som$name %in% colnames(frogs_in_ponds),]

# transform into presence-absence
species_for_som[species_for_som > 0] <- 1

# transform into long format
species_for_som_long <- melt(species_for_som, id.vars = "name", measure.vars = c(2:257))

# sample names into variables
pcr_replicate <- sapply(strsplit(as.character(species_for_som_long$variable), 
                                 split='.16S'), 
                        function(x) (x[2]))

pond_names_som <- substr(as.character(species_for_som_long$variable), 8, 9)

last_part_name <- substr(as.character(species_for_som_long$variable),
                         (nchar(as.character(species_for_som_long$variable))+2)-6,
                         nchar(as.character(species_for_som_long$variable)))

preservation_extraction <- substr(last_part_name,0,1)

# SOM long format with replicate descriptors
som_long_predictors <- data.frame(name = species_for_som_long$name,
                                  pond = pond_names_som,
                                  method = preservation_extraction,
                                  replicate = pcr_replicate)

write.csv(file="som_input_170622.csv", som_long_predictors)



# Results from Thierry
# SOM1: all species modelled together - for a general comparison of the 
# three preservation / extraction approaches.
# psi = Pr(occupancy)
# p = Pr(detection)
# fp = Pr(flase pos.)
# b = Pr(ambiguous detection occurs) == only 1 PCR positive / 4
SOM1 = read.csv(file = "SOM/SOM_results1.csv", row.names = 1, header=T)

# Plot method effects
# Detection probabilities per method
pDetect = grep('^p\\(meth.\\)', rownames(SOM1))
par(mfrow=c(2,1), mar=c(4,10,1,1))
plot(SOM1[pDetect,"estimate"], c(3:1), xlim=c(min(SOM1$X2.5..CI[pDetect]),
                                              max(SOM1$X97.5..CI[pDetect])), 
     yaxt = "n", xlab = "Detection probability", ylab="", 
     ylim=c(0.5,3.5), pch=19)
segments(SOM1[pDetect,2], c(3:1), SOM1[pDetect,3], c(3:1))
axis(2, at=c(3:1), label=c("GFF (2 um, in liquid: A)", 
                           "GFF (2 um, dried: B)",
                           "Nylon (0.2 um, dried: C)"), las=1)

# Plot method false positives
# False positives per method
pFPos = grep('^fp\\(meth.\\)', rownames(SOM1))
plot(SOM1[pFPos,"estimate"], c(3:1), xlim=c(min(SOM1$X2.5..CI[pFPos]),
                                              max(SOM1$X97.5..CI[pFPos])), 
     yaxt = "n", xlab = "False positives", ylab="", 
     ylim=c(0.5,3.5), pch=19)
segments(SOM1[pFPos,2], c(3:1), SOM1[pFPos,3], c(3:1))
axis(2, at=c(3:1), label=c("GFF (2 um, in liquid: A)", 
                           "GFF (2 um, dried: B)",
                           "Nylon (0.2 um, dried: C)"), las=1)

# SOM2: species modelled separately - for species-specific methodological differences
# Detection probabilities for each species per method
SOM2.p = read.csv(file = "SOM/SOM2_detection_prob.csv", row.names = 1, header=T)

# False positives, each species per method
SOM2.fp = SOM2.p = read.csv(file = "SOM/SOM2_false_pos.csv", row.names = 1, header=T)

# Coefficient plots
par(mfrow=c(1,3), mar=c(4.1,1,0.5,0.5))
plot(rep(0,nrow(SOM2.p)), c(nrow(SOM2.p):1), type="n", 
     bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
plot(SOM2.p[,"pA_estimate"], c(nrow(SOM2.p):1)-0.2, 
     xlim=c(0,0.5), yaxt="n",
     xlab="Detection probability", ylab="", pch=19, cex=0.5, 
     ylim=c(nrow(SOM2.p),1))
points(SOM2.p[,"pB_estimate"], c(nrow(SOM2.p):1), pch=17, cex=0.5)
points(SOM2.p[,"pC_estimate"], c(nrow(SOM2.p):1)+0.2, pch=15, cex=0.5)
segments(SOM2.p$pA_2.5..CI, c(nrow(SOM2.p):1)-0.2, 
         SOM2.p$pA_97.5..CI, c(nrow(SOM2.p):1)-0.2)
segments(SOM2.p$pB_2.5..CI, c(nrow(SOM2.p):1), 
         SOM2.p$pB_97.5..CI, c(nrow(SOM2.p):1))
segments(SOM2.p$pC_2.5..CI, c(nrow(SOM2.p):1)+0.2, 
         SOM2.p$pC_97.5..CI, c(nrow(SOM2.p):1)+0.2)
for(i in 1:nrow(SOM2.p)){
  lines(c(-0,1), c(i,i)-0.5, lty="dotted")
}
axis(2, at=c(nrow(SOM2.p):1), label=rownames(SOM2.p), las=1, cex.axis=0.9,
     tick=F)
plot(SOM2.fp[,"pA_estimate"], c(nrow(SOM2.p):1)-0.2, 
     xlim=c(0,0.5), yaxt="n",
     xlab="False positives", ylab="", pch=19, cex=0.5, 
     ylim=c(nrow(SOM2.p),1))
points(SOM2.fp[,"pB_estimate"], c(nrow(SOM2.fp):1), pch=17, cex=0.5)
points(SOM2.fp[,"pC_estimate"], c(nrow(SOM2.fp):1)+0.2, pch=15, cex=0.5)
segments(SOM2.fp$pA_2.5..CI, c(nrow(SOM2.fp):1)-0.2, 
         SOM2.fp$pA_97.5..CI, c(nrow(SOM2.fp):1)-0.2)
segments(SOM2.fp$pB_2.5..CI, c(nrow(SOM2.fp):1), 
         SOM2.fp$pB_97.5..CI, c(nrow(SOM2.fp):1))
segments(SOM2.fp$pC_2.5..CI, c(nrow(SOM2.fp):1)+0.2, 
         SOM2.fp$pC_97.5..CI, c(nrow(SOM2.fp):1)+0.2)
for(i in 1:nrow(SOM2.p)){
  lines(c(-0,1), c(i,i)-0.5, lty="dotted")
}





# Amphibians
Amphi = c("Dendropsophus leucophyllatus","Dendropsophus melanargyreus",
          "Dendropsophus minutus","Dendropsophus nanus","Dermatonotus muelleri",
          "Elachistocleis sp.","Eupemphix nattereri","Hylinae","Hyloidea",
          "Hypsiboas punctatus","Hypsiboas raniceps","Leptodactylus",
          "Leptodactylus chaquensis","Leptodactylus elenae","Leptodactylus fuscus",
          "Leptodactylus podicipinus","Leptodactylus syphax","Leptodactylus vastus",
          "Osteocephalus","Osteocephalus taurinus","Phyllomedusa azurea",
          "Phyllomedusa boliviana","Physalaemus albonotatus","Physalaemus centralis",
          "Physalaemus cuvieri","Pseudis limellum","Pseudis paradoxa",
          "Pseudopaludicola mystacalis","Rhinella schneideri","Scinax fuscomarginatus",
          "Scinax fuscovarius","Scinax nasicus","Scinax ruber","Scinax sp. FB-2014a",
          "Bufonidae","Elachistocleis","Gastrophryninae","Hylidae",
          "Leptodactylus latinasus","Melanophryniscus","Pseudis","Pseudis laevis","Scinax")
nrow(AbundSOM[Amphi,])

InputSOM = AbundSOM[Amphi,grep("sample.T", names(AbundSOM))]

# create SOM alternative input for export
cbind(rownames(InputSOM)[1], t(InputSOM[1,]))

ArrangedSOM = matrix(ncol = 2)
for (i in 1:nrow(InputSOM)) {
  ArrangedSOM = rbind(ArrangedSOM, cbind(rownames(InputSOM)[i], t(InputSOM[i,])))
}

# This writes out the file that can be edited to get the input format required by Thierry
# write.csv(file="input_som.csv", ArrangedSOM)

# Clarify vague identifications
# Osteocephalus
grep("Osteocephalus", MetaHead$sci_name)
osteo = MetaHead[grep("Osteocephalus", MetaHead$sci_name),]

# Elachistocleis
elach = MetaHead[grep("Elachistocleis", MetaHead$sci_name),]

# Pseudis
pseud = MetaHead[grep("Pseudis", MetaHead$sci_name),]

write.csv(file = "species-to-clear.csv", rbind(osteo,elach,pseud))

# tree checked in seaview, I would believe all variants.

# 
  
  
  
  
  
  
  








## Hominidae sequences
sum(apply(AbundHead, 1, sum)[MetaHead$family == "Hominidae"], na.rm=T)


summary(MetaAll$sci_name)


# Read distributions
hist(apply(AllData[,8:679], 2, sum), nclass = 20)
sum(AllData[,8:679])

# Calculate how often a sequence is head, internal or singleton
AllData = cbind(AllData, 
                no_sing = apply(AllData[,680:1351], 1, function(x) length(which(x == "s"))),
                no_inte = apply(AllData[,680:1351], 1, function(x) length(which(x == "i"))),
                no_head = apply(AllData[,680:1351], 1, function(x) length(which(x == "h"))))

HasHead = AllData$no_head != 0

# Keep sequence variants that were at least once observed as head 
DataWitHead = AllData[HasHead,]
sum(DataWitHead[,8:679])
hist(apply(DataWitHead[,8:679], 2, sum), nclass = 20)
summary(apply(DataWitHead[,8:679], 2, sum))

# retain only 16S samples
All_16S = grep("16S.$", names(DataWitHead))
Heads16S = cbind(DataWitHead[,1:7], DataWitHead[,All_16S], DataWitHead[,1352:1362])
sum(Heads16S[,8:343])

write.csv(file="16S_R_filtered.csv", Heads16S)





TotPresent = apply(LakeReads,2,function(vec) sum(vec>0))
hist(TotPresent)


###

# Map
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggsn)
library(rgdal)
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



