################################################################################
# this code carriers out analysis for type IIS endonucleases restriction sites
# in the coronavirus genome 
################################################################################

#### Load the needed libraries
library(Biostrings)
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(DECIPHER)
library(dplyr)
library(ggbeeswarm)
library(BuenColors)

#### set the data directory and import data
setwd("/Users/Downloads/")

# import all betacoronavirus genomes from Jan 1, 1900 to Nov 30, 2019, total 1316 complete genomes when the study began
# Similar analyses for alphacoronavirus and SARS-CoV-2 were based on this script with the corresponding genome data.
# Alphacoronavirus and betacoronavirus genome (complete) were downloaded from NCBI virus
# SARS-CoV-2 complete genomes were downloaded from GISAID (https://gisaid.org/).


vref = readDNAStringSet("betacornavirus_11.30.2019_1316.fasta") 

fastaFile <- vref
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)


### remove genomes without the BsmBI / BsaI sites
data(RESTRICTION_ENZYMES) # load the restriction enzymes
numb_cuts <- data.frame()

for (i in 1:length(df$sequence)) {
  cuts <- DigestDNA(RESTRICTION_ENZYMES[c("BsaI", "BsmBI")], df$sequence[i], type="positions", strand="top")
  temp_cuts <- length(unlist(cuts)) # how many cutting sites in total
  temp_cuts <- as.data.frame(temp_cuts)
  temp_cuts$seq_name <- df$seq_name[i]
  numb_cuts <- rbind(numb_cuts, temp_cuts)
}

BB.cuts <- numb_cuts[which(numb_cuts$temp_cuts > 0), ] 

df.3bb <- merge(df, BB.cuts, by = "seq_name")
df.3bb$noChar <- nchar(df.3bb$sequence)

# write.csv(df.3bb, "betacoronavirus_genome12.30.19_BsmBI_BsaI>2sites_1324.csv")


# find the recognition sites and cutting band size
data(RESTRICTION_ENZYMES) # load the restriction enzymes

re.site <-data.frame()

for (i in  1:length(df.3bb$sequence)) { 

# find the recognition sites
  cuts <- DigestDNA(RESTRICTION_ENZYMES[c("BsaI", "BsmBI")], df.3bb$sequence[i], type="positions", strand="top")
  length(unlist(cuts)) # how many cuts
  temp_sites <- as.data.frame(cuts[[1]])

## add the band size difference (keep the first top value, and use the second site - the first site)
  temp_sites$seq_name <- df.3bb$seq_name[i]
  temp_sites$noChar <- df.3bb$noChar[i] 
  temp_sites$Difference <- temp_sites$top - lag(temp_sites$top, default = 0)  

## add the last fragment, i.e the band after the last recognition enzyme site
  lastrow <- tail(temp_sites, n =1)
  lastrow$Difference <- lastrow$noChar - lastrow$top  # size is equal to the genome size minus the last RE site, + 1?
  lastrow$top[1] <-  NA # assign the lastrow top to NA
  temp_sites <- rbind(temp_sites, lastrow)
  
  re.site <- rbind(re.site, temp_sites)  

}  


###### merge with the original data
re.site2 <-  re.site  
re.site.fs.nb.mg.all <- merge(re.site2, df, by = "seq_name")


#### find the relative positions of BsmBI/BsaI on the genome
re.site.fs.nb.mg.all$rel.top <- re.site.fs.nb.mg.all$top / re.site.fs.nb.mg.all$noChar

# short the name of accession id
re.site.fs.nb.mg.all$Accessionid <- substr(re.site.fs.nb.mg.all$seq_name, 1, 9) 

# plot all the BsmBI/BsaI sites in the coronovirus genomes
ggplot(data =re.site.fs.nb.mg.all) +
  aes(y = rel.top, x = Accessionid) +
  geom_beeswarm(size = 0.5, alpha = 0.2) +
  coord_flip() +
  ylim(0,1) + theme_light()



#### find genomes met the IVGA criteria, i.e. 4-7 RE sites, and each <8000bp,
re.site.fs.nb <- re.site.fs.nb.mg.all %>% group_by(seq_name) %>% tally()

re.site.fs.nb2 <- re.site.fs.nb[which(re.site.fs.nb$n > 4 & re.site.fs.nb$n < 9), ] # n is the total number of fragments after BsmBI/BsaI digestion
re.site.fs.nb.mg <- merge(re.site.fs.nb2,re.site.fs.nb.mg.all, by = "seq_name")

# remove genomes with >8000 (8500) bp fragments after BsmBI/BsaI digestion
re.site.fs.nb.mg22 <- re.site.fs.nb.mg %>%
  group_by(seq_name) %>%
  filter(!any(Difference > 8000 ))  # 8000 or 8500 


#### plot the distance bw restriction sites across the genome level 
re.site.fs.nb.mg22$rel.top <- re.site.fs.nb.mg22$top /re.site.fs.nb.mg22$noChar

re.site.fs.nb.mg22$Accessionid <- substr(re.site.fs.nb.mg22$seq_name, 1, 9) 

ggplot(data = re.site.fs.nb.mg22) +
  aes(y = rel.top, x = Accessionid) +
  geom_beeswarm(color = "cornflowerblue") +
  coord_flip() +
  ylim(0,1) + theme_light()



#### hightlight those genomes
re.site.fs.nb.mg.all$label <- "other"
re.site.fs.nb.mg22$label <- re.site.fs.nb.mg22$Accessionid
re.site.fs.nb.mg22 <- re.site.fs.nb.mg22[, -2] 


#### Add sars-cov-2 
posctl <- read.csv("NC_045512_cutting sites.csv", header = TRUE)
posctl <- posctl[, -1]


# seperate the highlighted genomes and others, and merge with sars-cov-2 
re.site.fs.nb.mg.all.others <- re.site.fs.nb.mg.all[!(re.site.fs.nb.mg.all$seq_name %in% re.site.fs.nb.mg22$seq_name), ]

data.all1 <- rbind(re.site.fs.nb.mg22, posctl)

# plot
ggplot() +
  geom_point(data= re.site.fs.nb.mg.all.others, aes(y = rel.top, x = Accessionid), color = "grey50", size = 0.5, alpha = 0.1) +
  geom_vline(xintercept = c(unique(data.all1$label)), size = 0.2, color = "grey80")+
  geom_point(data = shuf(data.all1), aes(y = rel.top, x = Accessionid, color = label), size = 2) +
  coord_flip() +
  ylim(0,1) + theme_minimal()  + theme(panel.grid.major.y = element_line(colour="grey97", size = 0.05), 
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank())


#### for visualization purpose, manually change the relative position of those highlighted genomes in the plot
#### by renaming them based on the order in the matrix. 

c <- as.data.frame((re.site.fs.nb.mg.all.others$Accessionid))
data.all1.rename <- data.all1

# if using 5~8 fragments, < 8000 bp, then use the following to plot
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419105."] <- "JF292914.z" # c[1200, 1]
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419106."] <- "KF600632.z" # c[2600, 1]
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419107."] <- "KF923892.z" # c[2800, 1]
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419109."] <- "KF923913.z" # c[3000, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419111."] <- "KJ156934.z" # c[3200, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419112."] <- "KC776174.z" # c[1600, 1]
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419110."] <- "JX163923.z" # c[1400, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "KY419113."] <- "KJ650098.z" # c[3400, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MG762674."] <- "KY417144.z" # c[5600, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MN514964."] <- "MF598682.z" # c[7800, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MN514965."] <- "MF598704.z" # c[8000, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MN514967."] <- "MG011341.z" # c[8200, 1] 
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MW165134."] <- "MG977452.z" # c[9200, 1]
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MW202339."] <- "MN723543.z" # c[10500, 1]
data.all1.rename$Accessionid[data.all1.rename$Accessionid == "NC_045512"] <- "NC_045512" # c[11000, 1]

# # if using 5~8 fragments, <8500 bp, then add two more
# data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MW218395."] <- "HM211101.z" # c[762, 1]
# data.all1.rename$Accessionid[data.all1.rename$Accessionid == "MH002339."] <- "FJ882946.z" # c[400, 1] 

# plot
ggplot() +
  geom_point(data= re.site.fs.nb.mg.all.others, aes(y = rel.top, x = Accessionid), color = "grey50", size = 1, alpha = 0.03) +
  geom_vline(xintercept = c(unique(data.all1.rename$Accessionid)), size = 0.2, color = "grey80")+
  geom_point(data = data.all1.rename, aes(y = rel.top, x = Accessionid, color = label), size = 2, alpha = 1) +
  coord_flip() +
  ylim(0,1) + theme_minimal() + theme(panel.grid.major.y = element_line(colour="grey97", size = 0.05),
                                      panel.grid.minor = element_blank(),
                                      panel.background = element_blank())


