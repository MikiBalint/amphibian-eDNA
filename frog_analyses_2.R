library(vegan)
# library(vegan3d)
library(bvenn)
library(knitr)
library(boral)
library(mvabund)
# library(geosphere)
library(car)
# library(ape) # installed with ctv, infos here: http://www.phytools.org/eqg/Exercise_3.2/

OwnAssign = read.csv(file="abundances/frogs_16S_own.tab", sep="\t",
                     header=T, row.names=1)
EmblAssign = read.csv(file="abundances/frogs_16S_EMBL.tab", sep="\t",
                      header=T, row.names=1)

# Searches for particular taxa and species, example:
OwnAssign[OwnAssign$scientific_name == "Leptodactylus chaquensis",]

# Remove sequences from the 12S experiment. This expereriment 
# was multiplexed into the sequencing run of the current
# project, but it is not analyzed here.
EmblAssign = EmblAssign[EmblAssign[,"experiment"] == "16S",]

# Metadata
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
  
# remove the 12S statuses
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

# remove the 12S abundances
AbundEmbl = (AbundEmbl[, grepl("16", colnames(AbundEmbl))])
AbundAll = rbind(AbundOwn, AbundEmbl)

# Sequnce variant abundances at least once head
AbundHead = AbundAll[MetaAll$h_count > 0,]
MetaHead = MetaAll[MetaAll$h_count > 0,]

write.csv(file="AbundHead.csv", AbundHead)

# 160525 clean up negative controls: remove the maximum read number of a 
# sequence variant found in a negative control from every sample that 
# contains that sequence variant

# negative controls

# PCR controls
PNC = grep("sample.P.NC", names(AbundHead))
colnames(AbundHead)[PNC]

# glass fiber filter extraction control
NTC = grep("sample.NTC", names(AbundHead))
colnames(AbundHead)[NTC]

# nylon filter extraction control
NC = grep("sample.NC", names(AbundHead))
colnames(AbundHead)[NC]

# multiplexing control
MPX = grep("sample.MPX", names(AbundHead))
colnames(AbundHead)[MPX]

# Maximum number of reads in any control sample
MaxControl = apply(AbundHead[,c(PNC,NTC,NC,MPX)], 1, max)

# Extract the highest read number in a control from every sample
# checks
# AbundHead[1:5,100:102]
# MaxControl[1:5]
# sweep(AbundHead[1:5,100:102], 1, MaxControl[1:5], "-")
AbundControlled = sweep(AbundHead, 1, MaxControl, "-")

AbundControlled[AbundControlled < 0] <- 0

summary(apply(AbundControlled,2,sum))
summary(apply(AbundHead,1,sum))

# combine the replicates of samples
# get sample names, code from here: http://stackoverflow.com/questions/9704213/r-remove-part-of-string
SampleNames = levels(as.factor(sapply(strsplit(names(AbundControlled), 
                                               split='16S', fixed=TRUE), 
                                      function(x) (x[1]))))

# Sum the replicates for each sample
SummedReps = data.frame(row.names = rownames(AbundControlled))
for (i in 1:length(SampleNames)){
  ActualSet = grep(SampleNames[i], names(AbundControlled)) # grep the columns of interest
  SummedReps = cbind(SummedReps, apply(AbundControlled[ActualSet], 1, sum))
}
colnames(SummedReps) = SampleNames
# write.csv(file="test.csv", AbundControlled)

# In how many replicates observed per sample?
PresentReps = data.frame(row.names = rownames(AbundControlled))
for (i in 1:length(SampleNames)){
  ActualSet = grep(SampleNames[i], names(AbundControlled))
  Selected = AbundControlled[ActualSet]
  Selected[Selected > 0] <- 1 # set the read numbers to 1
  PresentReps = cbind(PresentReps, apply(Selected, 1, sum))
}
colnames(PresentReps) = SampleNames

# Set read numbers to 0 in a sample if the sequence variant was not observed in at least
# two PCR replicates
SummedControlled = SummedReps
SummedControlled[PresentReps < 2] <- 0

# Categories for sequence variants: frogs, humans, other aquatic stuff, etc.
# Amphi, HighGroup, "Homo sapiens", FarmAnim, Bird, Fish, Insect, Mammal
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

HighGroup = c("root", "Bilateria", "Amniota")

FarmAnim = c("Canis", "Sus scrofa", "Gallus", "Mus", "Bos", "Capreolus capreolus")
# write.csv(file = "other_taxa.csv", 
#           MetaHead$sci_name[!MetaHead$sci_name %in% 
#                               c(Amphi, HighGroup, "Homo sapiens", FarmAnim)])

Bird = c("Meleagris gallopavo", "Jacana", "Gallus")

Fish = c("Astyanax","Brachyplatystoma rousseauxii","Callichthys callichthys","Characidae",
         "Characiphysae","Opisthonema oglinum","Salmonidae","Sardinops","Sardinops sagax",
         "Siluriformes","Siluroidei","Synbranchus marmoratus","Thraupidae","Thunnini")

Insect = c("Micrathyria ocellata", "Libellulidae", "Micrathyria", "Tabanus", "Tramea")

Mammal = c("Homo sapiens", "Canis", "Sus scrofa", "Mus", "Bos", "Capreolus capreolus",
           "Myotis", "Laurasiatheria")

# Proportions of different groups
par(mar=c(1,1,1,1))
pie(c(sum(SummedControlled[MetaHead$sci_name %in% Amphi,]), 
      sum(SummedControlled[MetaHead$sci_name %in% Mammal,]),
      sum(SummedControlled[MetaHead$sci_name %in% Fish,]),
      sum(SummedControlled[MetaHead$sci_name %in% Insect,]),
      sum(SummedControlled[MetaHead$sci_name %in% Bird,]),
      sum(SummedControlled[MetaHead$sci_name %in% HighGroup,])),
    labels = paste(c("Frogs\n2 198 212", "Mammals\n38 568", 
                     "Fish\n1 873 599", "Insects\n314 910", 
                     "Birds\n4 785", "Higher groups\n1 384 940")),
    col = c(gray(0.9), gray(0.75), gray(0.60), gray(0.45), gray(0.30), gray(0.15)))

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
FrogAggregate = data.frame(name = MetaHead$sci_name[MetaHead$sci_name %in% Amphi], 
                         SummedControlled[MetaHead$sci_name %in% Amphi,])
FrogCounts = aggregate(. ~ name, FrogAggregate, sum, na.action = na.exclude)
rownames(FrogCounts) = FrogCounts$name
FrogCounts = FrogCounts[,2:96]

# Correct the Scinax madeirae: Scinax sp. FB-2014a is actually S. madeirae
# Species list
rownames(FrogCounts)[34] <- "Scinax madeirae"

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
Lakes = substring(names(FrogCounts), 8, 9)

# Site metadata: lakes, preservation/extraction methods, read numbers
SiteMeta = data.frame(Lakes = Lakes, 
                      Methods = c(rep(c("Control"),31), 
                                  rep(c("A","B","C"),12),
                                  rep("A",4), rep("B",4), rep("C",3),
                                  rep("A",3), rep("B",3), rep("C",3),
                                  rep("A",3), rep("B",3), rep("C",2)),
                      Reads = apply(FrogCounts, 2, sum))
LakeCodes = c("T1", "T2", "T3", "T4", "T5")

# Grep names of all samples. Sample names are in the 'lakeCodes'
IsLakeSample = grep(paste(LakeCodes, collapse='|'), 
                    names(FrogCounts), ignore.case=TRUE)

# Create read abundance table for frog species found in the ponds
FrogsInLakes = FrogCounts[,IsLakeSample]
FrogsInLakes = FrogsInLakes[apply(FrogsInLakes,1,sum) > 0,]

# Sum up species counts by lakes
FrogLakes = data.frame(row.names = rownames(FrogsInLakes))
for (i in 1:length(LakeCodes)) {
  ActualSet = grep(LakeCodes[i], names(FrogsInLakes))
  FrogLakes = cbind(FrogLakes, apply(FrogsInLakes[ActualSet], 1, sum))
}
colnames(FrogLakes) = LakeCodes

write.csv(file="Frog_sums_ponds.csv", 
          cbind(FrogLakes, Total = apply(FrogLakes,1,sum)))

# Total eDNA detected species in a lake
apply(FrogLakes,2,function(x) sum(x>0))

# Get positive samples + the field negative tapwater sample
FieldCodes = c("Cen", "Lvastus")
PositiveCodes = c("PCE", "PCUNE")

# Grep names of positive and field tap water samples. 
# Sample names are in the 'lakeCodes'
IsPosField = grep(paste(c(FieldCodes, PositiveCodes), collapse='|'), 
                    names(FrogCounts), ignore.case=TRUE)

# Create read abundance table for frog species found in the ponds
FrogsInPosField = FrogCounts[,IsPosField]
FrogsInPosField = FrogsInPosField[apply(FrogsInPosField,1,sum) > 0,]

# Sum up species counts by positive and field controls
FrogPosField = data.frame(row.names = rownames(FrogsInPosField))
for (i in 1:length(c(FieldCodes, PositiveCodes))) {
  ActualSet = grep(c(FieldCodes, PositiveCodes)[i], names(FrogsInPosField))
  FrogPosField = cbind(FrogPosField, apply(FrogsInPosField[ActualSet], 1, sum))
}
colnames(FrogPosField) = c(FieldCodes, PositiveCodes)

# 
write.csv(file="Frog_sums_FieldPos.csv", FrogPosField)

# Total eDNA detected species in a lake
apply(FrogPosField,2,function(x)sum(x > 0))

#######
# Venn diagrams (http://stackoverflow.com/questions/11722497/how-to-make-venn-diagrams-in-r)
SurveyList = c("Dendropsophus leucophyllatus",
               "Dendropsophus melanargyreus",
               "Dendropsophus minutus",
               "Dendropsophus nanus",
               "Dendropsophus salli",
               "Elachistocleis sp.",
               "Eupemphix nattereri",
               "Hypsiboas geographicus",
               "Hypsiboas punctatus",
               "Hypsiboas raniceps",
               "Leptodactylus fuscus",
               "Leptodactylus syphax",
               "Leptodactylus vastus",
               "Leptodactylus chaquensis",
               "Osteocephalus taurinus",
               "Phyllomedusa azurea",
               "Phyllomedusa boliviana",
               "Physalaemus albonotatus",
               "Physalaemus centralis",
               "Pseudis paradoxa",
               "Pseudis laevis",
               "Pseudopaludicola mystacalis",
               "Sphaenorhynchus lacteus",
               "Scinax fuscomarginatus",
               "Scinax fuscovarius",
               "Scinax madeirae",
               "Scinax ruber",
               "Scinax nasicus")

InWaterList = read.csv(file = "abundances/frogs_in_water.csv", header = F)
InWaterList = levels(c(InWaterList)$V1)

eDNAList = rownames(FrogLakes)[apply(FrogLakes,1,sum) > 0]

# eDNA and audio-visual survey list
bvenn(list(eDNA = eDNAList, "Audio-visual\nsurvey" = SurveyList), 
      fontsize=10)
# eDNA and in water frog list
bvenn(list(eDNA = eDNAList, "Frogs or tadpoles\nin water" = InWaterList),
      fontsize = 10)

# Positive controls
# PCE: DNA from 12 species in equal concentrations
C.PCE = grep("sample.P.PCE", names(FrogCounts))

# PCUNE: DNA from 12 species in stepwise doubled concentrations
C.PCUNE = grep("sample.P.PCUNE", names(FrogCounts))

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
                     "Dendropsophus nanus"))

# One of the 12 species (Hypsiboas geographicus) is not present in the PCE results.
IsInPCE = cbind(PosList, PosList %in% rownames(FrogCounts))

# Visualize conncentrations in the PCE. Lines are PCE samples.
palette(colors())
par(mar=c(8,4,1,1))
plot(c(1:12), seq(0.1, (max(FrogCounts[,C.PCE])), 
                  (max(FrogCounts[,C.PCE]))/12), type="n",
     xaxt="n", xlab="", ylab="Read numbers") 
axis(1, at = c(1:12), labels = PosList, las = 2, cex.axis=0.6)
for (i in C.PCE) {
  lines(c(1:12), (FrogCounts[PosList,i]), col=i+10, lwd=2)
}

# The seven PCE samples are positively correlated at R > 0.7.

# Correction factors for PCE (obsolate)
PCEConc = 5 #ng/ul
PCEUNEConc = read.csv(file="PCUNE_conc.csv", header=T, row.names = 1)
# 
# #  the mean abundance of the species in PCE with 1x dilution 
# # in PCUNE should correspond the  5 ng/ul conc: Phyllomedusa azurea
# mean(as.numeric(FrogCounts["Phyllomedusa azurea", C.PCE]))
# 
# # Divide this by all abundances of control species to get the correction factor
# CorrectFactor = mean(as.numeric(FrogCounts["Phyllomedusa azurea", C.PCE])) /
#   apply(FrogCounts[PosList,],1,mean)
# 
# # Corrected PCE abundances
# FrogCountsCorrected = FrogCounts[PosList,] * CorrectFactor
# apply(FrogCountsCorrected, 1, sum)

# Evaluation of non-equimolar concentrations and read numbers
# Correlations with the original DNA concentrations
# PCUNE abundances VS original concentrations

RangeX = range(PCEUNEConc[,"conc"])
RangeY = range(apply(FrogCounts[,C.PCUNE],1,max))

par(mfrow = c(1,1), mar=c(4,4,3,1))
plot(RangeX, RangeY, type="n", xlab="DNA concentrations", ylab = "Read numbers", 
     log="x", main = "Positiv non-equimolar controls")
for (i in C.PCUNE){
  points(jitter(PCEUNEConc[PosList[PosList != "Hypsiboas geographicus"],"conc"]), 
       FrogCounts[PosList[PosList != "Hypsiboas geographicus"],i], 
       lwd=2, col = i)
  # abline(lm(FrogCounts[PosList[PosList != "Hypsiboas geographicus"],i] ~ 
  #             PCEUNEConc[PosList[PosList != "Hypsiboas geographicus"],"conc"]), 
  #        col = i)
}

#PCUNE correlations
PCUNE.cor = c()
for (i in C.PCUNE){
  PCUNE.cor = c(PCUNE.cor, cor(PCEUNEConc[PosList[PosList != "Hypsiboas geographicus"],"conc"], 
         FrogCounts[PosList[PosList != "Hypsiboas geographicus"],i]))
}
mean(PCUNE.cor)

# Correlation of concentrations and read numbers in the PCUNE
ConcMeanRead = data.frame()
for (i in PosList[PosList != "Hypsiboas geographicus"]){
  ConcRead = c(PCEUNEConc[i,"conc"], mean(as.numeric(FrogCounts[i,C.PCUNE])),
               mean(as.numeric(FrogCounts[i,C.PCE])))
  ConcMeanRead = rbind(ConcMeanRead, ConcRead)
}
colnames(ConcMeanRead) = c("conc", "MeanReadPCUNE", "MeanReadPCE")
rownames(ConcMeanRead) = PosList[PosList != "Hypsiboas geographicus"]

cor(ConcMeanRead, method = "pearson")


# # Visualize conncentrations in the PCE. Lines are PCE samples.
# palette(colors())
# par(mar=c(8,4,1,1))
# plot(c(1:11), seq(0.1, (max(FrogCountsCorrected[,C.PCE])), 
#                   (max(FrogCountsCorrected[,C.PCE]))/11), type="n",
#      xaxt="n", xlab="", ylab="Read numbers") 
# axis(1, at = c(1:11), labels = PosList[PosList != "Hypsiboas geographicus"], las = 2, cex.axis=0.6)
# for (i in C.PCE) {
#   lines(c(1:11), (FrogCountsCorrected[PosList != "Hypsiboas geographicus",i]), col=i+10, lwd=2)
# }



#####
# Ecological signal: differences among the lakes VS preservation/extraction methods
# Transpose the frog counts: species in columns, samples in lines
FrogCountsT = t(FrogCounts)

# Filter for lake samples
Samples = SiteMeta$Methods != "Control"
FrogCountsLake = FrogCountsT[Samples,]


# Some samples don't have reads at all, some species were not seen in lakes
summary(apply(FrogCountsLake,1,sum))
summary(apply(FrogCountsLake,2,sum))

# Meaningful species matrix for lakes
FrogCountsNoZero = FrogCountsLake[apply(FrogCountsLake,1,sum) > 0,
                                  apply(FrogCountsLake,2,sum) > 0]

# Corresponding metadata
MetaLakesNoZero = SiteMeta[rownames(SiteMeta) %in% rownames(FrogCountsNoZero),]

# Remove species present in only one sample
FrogCountsNoOne = FrogCountsNoZero[,apply(FrogCountsNoZero,2,function(vec) sum(vec>0)) > 1]

# still species in every sample?
summary(apply(FrogCountsNoOne,1,sum))

# Model-based comparison of the lakes
FrogsMvabund = mvabund(FrogCountsNoOne)

# Multispecies GLM and ANOVA. Reads control for differences in sequencing depth.
m1 = manyglm(FrogsMvabund ~ Reads + Lakes, 
             data=MetaLakesNoZero)

# Increase nBoot=1000, and include p.uni for species-level statistics
m1.anova = anova(m1, nBoot = 10, p.uni = "adjusted")

# nice anova table
kable(m1.anova$table)

# Visualize the community results
# Model-based ordination with reads as covariates

# Overdispersion parameter
summary(m1$theta)
hist(m1$theta)

set.prior = set.prior <- list(type = c("normal","normal","normal","uniform"),
                              hypparams = c(100, 20, 100, 1))

# Unconstrained ordination
# ord.m.noconst <- boral(FrogCountsNoOne, family = "negative.binomial", num.lv = 2,
#                        prior.control = set.prior,
#                        n.burnin = 10, n.iteration = 100, n.thin = 1)
# 
# par(mfrow = c(2,2))
# plot(ord.m.noconst)

# Constrained ordination with sequencing depth
# Just keep the old 2 LV plot
ord.m.const <- boral(FrogCountsNoOne, X=MetaLakesNoZero$Reads, 
                         prior.control = set.prior,
                         family = "negative.binomial", num.lv = 2, 
                         n.burnin = 10, n.iteration = 100, n.thin = 1)

# Diagnostics
par(mfrow = c(2,2))
plot(ord.m.const)

# Colors for ordination
palette(colors())
LakeCol = gsub("T1", "green", MetaLakesNoZero$Lakes)
LakeCol = gsub("T2", "orange", LakeCol)
LakeCol = gsub("T3", "red", LakeCol)
LakeCol = gsub("T4", "purple", LakeCol)
LakeCol = gsub("T5", "blue", LakeCol)
MyColors = data.frame(cols = c("green","orange","red","purple","blue"),
                      lakes = levels(factor(MetaLakesNoZero$Lakes)))

# Constrained ordination plot
par(mfrow=c(2,1), mar=c(1,5,1,1), oma=c(1,1,0,0))
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE),
       heights=c(2,1))
ordiplot(ord.m.const$lv.median, type = "none", cex =2, 
                   display = "sites", # xlim = c(-0.1,0.1), ylim = c(-6,7),
                   main="Constrained LV ordination", 
                   cex.lab = 1.5, cex.axis = 1.5)
points(ordifrog,"sites", pch=20 ,col=LakeCol, cex = 2)
legend(5, 0, paste(c("T1", "T2", "T3", "T4", "T5")), fill=as.vector(MyColors$cols), 
       border="white", bty="n", cex = 1.5)
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
                         draw="polygon", col="green",
                         alpha=100,kind="se",conf=0.95, 
                         show.groups=(c("T1")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="orange",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T2")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="red",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T3")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="purple",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T4")))
ordiellipse(ordifrog, factor(MetaLakesNoZero$Lakes), cex=.5, 
            draw="polygon", col="blue",
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("T5")))

# Cluster of ponds using species presence-absence based on the audio-visual survey
VAES = read.csv(file = "abundances/VAES_presence-absence.csv", 
                header = T, row.names = 1)

# plot(hclust(vegdist(t(VAES), method = "jaccard"), method = "complete"))
par(mar = c(1,5,1,1))
plot(hclust(vegdist(t(VAES), method = "raup"), method = "average"),
     xlab = "", main = "", cex.lab = 1.5, cex.axis = 1.5, lwd = 1.5, cex = 1.5)
# plot(hclust(vegdist(t(VAES), method = "raup"), method = "average"))

# Unconstrained ordination plot
# plot(ord.m.noconst$lv.median, col=LakeCol, 
#      pch=19, main="Unconstrained LV ordination", las=1, xlab = "LV1", ylab = "LV2")
# legend(-6.5, 2, c("Small Croco", "Vastus", "Wetland basin", "Lacha Susanna", 
#                   "Wetland Centro"), fill=as.vector(MyColors$cols), 
#        border="white", bty="n")
# for (i in levels(factor(MetaLakesNoZero$Lakes))) {
#   ellipse(center = c(mean(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,"LV1"]), 
#                      mean(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,"LV2"])),
#           shape = cov(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,]),
#           radius = 0.95*mean(std.error(ord.m.noconst$lv.median[MetaLakesNoZero$Lakes == i,])), 
#           add=T, col=as.vector(MyColors$cols[MyColors$lakes == i]))
# }


# Distribution of read abundances
hist(apply(FrogCountsT,2,sum), nclass=20, col="grey", 
     main="Distribution of read abundances",
     xlab="Read counts per species", ylab="Frequency")

# Simple NMDS plot
# LakeNMDS = metaMDS(FrogCountsT[SamplesNoZero & Samples,])
# 
# par(mar=c(4,4,1,1))
# plot(LakeNMDS$points, type="n", xlab="NMDS1", ylab="NMDS2")
# ordispider(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"], col="grey")
# points(LakeNMDS, pch=20)
# mylegend = legend(-3.3, 3, c("Small Croco", "Vastus", "Wetland basin", "Lacha Susanna", "Wetland Centro"), 
#                   fill=c("orange","green","purple","red","blue"), border="white", bty="n")
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", 
#             col=c("orange"), alpha=170,kind="se",conf=0.95, show.groups=(c("T1")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("green"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T2")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("purple"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T3")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("red"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T4")))
# ordiellipse(LakeNMDS, SiteMeta[SamplesNoZero & Samples,"Lakes"],cex=.5, draw="polygon", col=c("blue"),
#             alpha=170,kind="se",conf=0.95, show.groups=(c("T5")))
# 
# 
# 
# 



# Preservation-extraction method evaluation with species occupancy models
# Script from Thiery comes here.

# Format the data for SOM
# !!! Should be updated.
# data frame of non-control samples
# positive contols: 
# C.PCE = grep("sample.P.PCE", names(SpeCounts))
# C.PCUNE = grep("sample.P.PCUNE", names(SpeCounts))
# 
# # negative controls: 
# C.PNC = grep("sample.P.NC", names(SpeCounts))
# C.NTC = grep("sample.NTC", names(SpeCounts))
# C.NC = grep("sample.NC", names(SpeCounts))
# C.MPX = grep("sample.MPX", names(SpeCounts))
# 
# Controls = c(C.PCE,C.PCUNE,C.PNC,C.NTC,C.NC,C.MPX)
# 
# # reads from lake observations
# AbundLakes = SpeCounts[-Controls]
# 
# # Transform to presence-absences for species occupancy models
# AbundSOM = AbundLakes
# AbundSOM[AbundSOM > 0] <- 1
# 
# # write.csv(file = "SOM_data.csv", AbundSOM)

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


