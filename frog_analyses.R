library(vegan)
library(vegan3d)
library(bvenn)

OwnAssign = read.csv(file="abundances/frogs_16S_own.tab", sep="\t",
                     header=T, row.names=1)
EmblAssign = read.csv(file="abundances/frogs_16S_EMBL.tab", sep="\t",
                      header=T, row.names=1)

# Remove sequences from the 12S experiments
EmblAssign = EmblAssign[EmblAssign[,"experiment"] == "16S",]

# Metadata
MetaOwn = data.frame(best_id = OwnAssign$best_identity.16S_fasta_nonredundant_noprimers,
                best_match = OwnAssign$best_match, count = OwnAssign$count, family = OwnAssign$family_name,
                genus = OwnAssign$genus_name, sci_name = OwnAssign$scientific_name, 
                taxid = OwnAssign$taxid, sequence = OwnAssign$sequence,
                seq_names = rownames(OwnAssign), row.names=9)
MetaEmbl = data.frame(best_id = EmblAssign$best_identity.db_16S,
                      best_match = EmblAssign$best_match, count = EmblAssign$count, family = EmblAssign$family_name,
                      genus = EmblAssign$genus_name, sci_name = EmblAssign$scientific_name, 
                      taxid = EmblAssign$taxid, sequence = EmblAssign$sequence,
                      seq_names = rownames(EmblAssign), row.names=9)
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

# Abundances
AbundOwn = OwnAssign[, grepl("sample", colnames(OwnAssign))]
AbundEmbl = EmblAssign[, grepl("sample", names(EmblAssign))]

# remove the 12S abundances
AbundEmbl = (AbundEmbl[, grepl("16", colnames(AbundEmbl))])
AbundAll = rbind(AbundOwn, AbundEmbl)

# Sequnce variant abundances at least once head
AbundHead = AbundAll[MetaAll$h_count > 0,]
MetaHead = MetaAll[MetaAll$h_count > 0,]

# Presence-absence data for the occupancy models

# 160525
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
MaxControl = summary(apply(AbundHead[,c(PNC,NTC,NC,MPX)], 1, max))

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

######################################################
# Deal with multiplexing and negative controls here!!!
######################################################

# Aggregate species
ToAggregate = cbind(names = Heads16S$scientific_name, Heads16S[,8:343])
SpeCounts = aggregate(. ~ names, ToAggregate, sum)
rownames(SpeCounts) = SpeCounts$names
SpeCounts = SpeCounts[,2:337]
# Species list
SpecList = sort(unique(rownames(SpeCounts)))

# Correct the Scinax madeirae 
rownames(SpeCounts)[32] <- "Scinax madeirae"

# Simple lake-based ordination
IsLake = grep("^T", colnames(SpeCounts))

# Get and transpose the abundances
LakeReads = t(SpeCounts[,IsLake])
hist(apply(LakeReads,2,sum), nclass=20, col="grey", main="Fajabundancia-eloszlás",
     xlab="Egyes fajok szekvencia-észlelései", ylab="Gyakoriság")

LakeMeta = substr(rownames(LakeReads),1,2)

boxplot(apply(LakeReads,1,sum) ~ LakeMeta)

LakeNMDS = metaMDS(LakeReads, k=3, trymax=100)

par(mar=c(4,4,1,1))
plot(LakeNMDS$points, type="n", xlab="NMDS1", ylab="NMDS2")
ordispider(LakeNMDS, LakeMeta, col="grey")
points(LakeNMDS, pch=20, cex=log(apply(LakeReads,1,mean))/3, col="black")
mylegend = legend(-1.6, 1.2, c("Small Croco", "Vastus", "Wetland basin", "Lacha Susanna", "Wetland Centro"), 
                  fill=c("orange","green","purple","red","blue"), border="white", bty="n")
ordiellipse(LakeNMDS, LakeMeta,cex=.5, draw="polygon", col=c("orange"),
                          alpha=170,kind="se",conf=0.95, show.groups=(c("T1")))
ordiellipse(LakeNMDS, LakeMeta,cex=.5, draw="polygon", col=c("green"),
            alpha=170,kind="se",conf=0.95, show.groups=(c("T2")))
ordiellipse(LakeNMDS, LakeMeta,cex=.5, draw="polygon", col=c("purple"),
            alpha=170,kind="se",conf=0.95, show.groups=(c("T3")))
ordiellipse(LakeNMDS, LakeMeta,cex=.5, draw="polygon", col=c("red"),
            alpha=170,kind="se",conf=0.95, show.groups=(c("T4")))
ordiellipse(LakeNMDS, LakeMeta,cex=.5, draw="polygon", col=c("blue"),
            alpha=170,kind="se",conf=0.95, show.groups=(c("T5")))

# # Dynamic plot
# ordirgl(LakeNMDS)
# orglspider(LakeNMDS, LakeMeta)
# orgltext(MDS.all.3, rownames(fun.some))

# Positiv equimolar controls
IsPCE = grep("^..PCE", colnames(SpeCounts))
PCEReads = t(SpeCounts[,IsPCE])

# Species in the positiv control
ConSpecList = sort(c("Phyllomedusa azurea",
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

# Not all species from the control are present in the dataset!!!
IsInPCE = cbind(ConSpecList, ConSpecList %in% SpecList)

# Load colors
palette(colors())

par(mar=c(8,4,1,1))
plot(c(1:12), seq(0.1, log(max(SpeCounts[,IsPCE])), log(max(SpeCounts[,IsPCE]))/12), type="n",
     xaxt="n", xlab="", ylab="log(szekvenciák száma)", ylim=c(3,10))
axis(1, at = c(1:12), labels = ConSpecList, las = 2, cex.axis=0.6)
for (i in IsPCE) {
  lines(c(1:12), log(SpeCounts[ConSpecList,i]), col=i, lwd=2)
}

TotPresent = apply(LakeReads,2,function(vec) sum(vec>0))
hist(TotPresent)

# Venn diagrams (http://stackoverflow.com/questions/11722497/how-to-make-venn-diagrams-in-r)
RegioList = c("Ameerega picta",
#####           
           "Ceratophrys sp.",
           "Chiasmocleis albopunctata",
           "Dendropsophus cachimbo",
           "Dendropsophus leali",
           "Dendropsophus leucophyllatus",
           "Dendropsophus melanargyreus",
           "Dendropsophus minutus",
           "Dendropsophus nanus",
           "Dendropsophus salli",
           "Dermatonotus muelleri",
           "Elachistocleis sp.",
           "Eupemphix nattereri",
           "Hypsiboas geographicus",
           "Hypsiboas punctatus",
           "Hypsiboas raniceps",
           "Leptodactylus cf. didymus",
           "Leptodactylus cf. diptyx",
           "Leptodactylus elenae",
           "Leptodactylus fuscus",
           "Leptodactylus leptodactyloides",
           "Leptodactylus mystacinus",
           "Leptodactylus syphax",
           "Leptodactylus vastus",
           "Oreobates heterodactylus",
           "Osteocephalus taurinus",
           "Phyllomedusa azurea",
           "Phyllomedusa boliviana",
           "Physalaemus centralis",
           "Physalaemus albonotatus",
           "Physalaemus cuvieri",
           "Pseudis paradoxa",
           "Pseudopaludicola mystacalis",
           "Rhinella cf. paraguayensis",
           "Rhinella schneideri",
           "Rhinella mirandaribeiroi",
           "Scinax fuscomarginatus",
           "Scinax fuscovarius",
           "Scinax madeirae",
           "Scinax nasicus",
           "Scinax ruber",
           "Sphaenorhynchus lacteus",
           "Trachycephalus cf. typhonius")
#####           

SurveyList = c("Dendropsophus leucophyllatus",
#####           
               "Dendropsophus melanargyreus",
               "Dendropsophus minutus",
               "Dendropsophus nanus",
               "Dendropsophus salli",
               "Elachistocleis sp.",
               "Eupemphix nattereri",
               "Hypsiboas geographicus",
               "Hypsiboas punctatus",
               "Hypsiboas raniceps(d)",
               "Leptodactylus fuscus",
               "Leptodactylus syphax",
               "Leptodactylus vastus",
               "Phyllomedusa azurea",
               "Phyllomedusa boliviana",
               "Physalaemus albonotatus",
               "Physalaemus centralis",
               "Pseudis paradoxa(d)",
               "Pseudopaludicola mystacalis",
               "Sphaenorhynchus lacteus",
               "Scinax fuscomarginatus",
               "Scinax fuscovarius",
               "Scinax madeirae",
               "Scinax ruber",
               "Scinax madeirae",
               "Scinax nasicus")
#####           

# eDNA and regional species list
bvenn(list(eDNS = SpecList, "Regionális felmérés" = RegioList), fontsize=18)

# eDNA and survey list
bvenn(list(eDNS = SpecList, "Helyi felmérés" = SurveyList), fontsize=18)

