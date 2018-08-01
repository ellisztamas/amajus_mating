# Import the monster AR dataset
core <- read.csv('~/Documents/Antirrhinum/core_data/FinalGenoEcol_short_nuc.csv', stringsAsFactors = F)

# L0057 is not recorded as alive in 2012, but it must have been because I collected seed from it.
which(sapply(strsplit(core$RametIDs, ","), function(x) "L0057" %in% x))
core$AliveRec_2012[core$PlantID_final == "L0057"] <- 1

# GPS positions
# pull out plants from Planoles in 2012
as2012<- core[core$AliveRec_2012 == 1 & core$Population == "Planoles",]
# Remove plants with odd names.
#These were included in datasets for cline analysis. Genotyping and error rates will be different, so remove these.
as2012 <- as2012[!as2012$PlantID_final %in% c("KTPla1","KTPla11","KTPla12", "KTPla15", "KTPla16", "KTPla18", "KTPla2","KTPla5",  "KTPla6", "KTPla7",  "KTPla8"),]

# There are some plants that are quite far into  the red flanks. We probably don't need these.
plot(as2012$Easting, as2012$Northing)
as2012 <- as2012[as2012$Easting < 426000,]

gps <- as2012[,c("PlantID_final","Easting", "Northing", "Altitude")]
# Centre Easting and Northing
gps$Easting  <- gps$Easting - mean(gps$Easting)
gps$Northing <- gps$Northing - mean(gps$Northing)
write.csv(gps, 'data/GPS_positions.csv', row.names = F)

sum(as2012$phenoCat_final == -9) # there are 358 plants with no phenotype data
# 13 genets have constituent ramets with conflicting phenotypes.
sum(as2012$phenoCat_final %in% c("WO;Y", "FR;FO","WR;WO","FR;WR;FR", "WR;FO;WR", "WO;WR","FR;FO;FR"))
# Assign Rosea genptypes where possible.
as2012$rosea[as2012$phenoCat_final %in% c("FR", "FO", "FR;FO", "FR;FO;FR")] <- "R/R"
as2012$rosea[as2012$phenoCat_final %in% c("WR", "WO", "WR;WO", "WO;WR")] <- "R/r"
as2012$rosea[as2012$phenoCat_final %in% c("W", "Y")] <- "r/r"
sum(is.na(as2012$rosea)) # 361 plants not assigned.
# Assign Sulfurea genotypes where possible
as2012$sulf[as2012$phenoCat_final %in% c("FR", "WR", "W")] <- "S/+"
as2012$sulf[as2012$phenoCat_final %in% c("FO", "WO", "Y", "WO;Y")] <- "s/s"
sum(is.na(as2012$sulf)) # 370 plants not assigned.
# write to disk
ros_sulf <- data.frame(id = as2012$PlantID_final, Ros = as2012$rosea, Sulf=as2012$sulf)
write.csv(ros_sulf, 'data/rosea_sulfurea.csv')

# import offspring data. They are in two files from separate typing efforts
offs1 <- read.csv('data/offspring_SNP_raw_data_1.csv', stringsAsFactors = F)
offs2 <- read.csv('data/offspring_SNP_raw_data_2.csv', stringsAsFactors = F)
# IDs match
all(offs1$ID == offs2$ID)
# bind into a single data frame.
offs <- cbind(offs1, offs2[-1])
rm(offs1, offs2)
# convert to FAPS format
source('R/lgc2faps.R')
offs <- lgc2faps(offs)

# pull out the names of the mothers.
offs$mothers <- sapply(strsplit(offs$ID, "_"), function(x) x[[1]])
# remove two maternal familes for whom we don't have maternal DNA. L0803.7 is from the array population Lys.
offs <- offs[! offs$mothers %in% c("K2063", "M0065", "L0803.7"),]
# Line up maternal names with the unqiue names assigned in core data.
# find the positions of the maternal names in the list of Ramets
mx <- sapply(offs$mothers, function(y) which(sapply(strsplit(as2012$RametIDs, ","), function(x) y %in% x)))
offs$mothers <- as2012$PlantID_final[mx]

# line up SNP names between parents and offspring
snp_summary <- read.csv('data/shortSummarySNPsAdj.csv', stringsAsFactors = F)
colnames(offs) <- c("ID", na.omit(snp_summary$LocusName[match(colnames(offs), snp_summary$LocusNameOLD)]), "mothers")

px <- intersect(colnames(offs), colnames(as2012)) # common set of SNPs.
parents <- lgc2faps(as2012[,c("PlantID_final",px)], "/")
offs <- offs[,c("ID","mothers",px)]

# write genotypes to disk.
write.csv(parents, file="data/parents_2012_genotypes.csv", row.names = F)
write.csv(offs,    file="data/offspring_2012_genotypes.csv", row.names = F)
