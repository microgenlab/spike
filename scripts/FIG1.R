library(stringr)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(viridis)
library(ggsignif)

#setwd("/home/csalazar/Documentos/1.spike/revisiones_spike_feb23/final_version/upload_github/")

#######################################################################
##################STANDARD-S METHOD####################################
#######################################################################

df <- read.table("./inputs/stdS_table.tsv", sep = "\t", header = T)
head(df)
unique(df$method)
df$depth <- as.numeric(df$depth)
df$completeness <- as.numeric(df$completeness)


#Convert all sublineages to a lineages/lineage-set
unique(df$genome_lineage)
df <- df %>%
  mutate(lineage_genome = case_when(
    endsWith(genome_lineage, "P.2") ~ "P.2",
    endsWith(genome_lineage, "P.6") ~ "P.6",
    endsWith(genome_lineage, "C.37.1") ~ "C.37",
    endsWith(genome_lineage, "P.1") ~ "P.1.*",
    endsWith(genome_lineage, "B.1.617.2") ~ "B.1.617.2",
    endsWith(genome_lineage, "B.1.351") ~ "B.1.351",
    endsWith(genome_lineage, "B.1.1.7") ~ "B.1.1.7",
    endsWith(genome_lineage, "AY.26") ~ "AY.*",
    endsWith(genome_lineage, "AY.122") ~ "AY.*",
    endsWith(genome_lineage, "B.1.621") ~ "B.1.621",
    endsWith(genome_lineage, "B.1.621.1") ~ "B.1.621",
    endsWith(genome_lineage, "C.37") ~ "C.37",
    endsWith(genome_lineage, "AY.20") ~ "AY.*",
    endsWith(genome_lineage, "AY.25.1") ~ "AY.*",
    endsWith(genome_lineage, "AY.99.2") ~ "AY.*",
    endsWith(genome_lineage, "BA.1") ~ "BA.1.*",
    endsWith(genome_lineage, "BA.1.1") ~ "BA.1.*"
  ))

unique(df$set_description)
df <- df %>%
  mutate(lineage_set = case_when(
    endsWith(set_description, "P.2") ~ "P.2",
    endsWith(set_description, "P.6") ~ "P.6",
    endsWith(set_description, "C.37.1") ~ "C.37",
    endsWith(set_description, "P.1_1") ~ "P.1.*",
    endsWith(set_description, "AY.30") ~ "AY.*",
    endsWith(set_description, "A_23") ~ "A_23",
    endsWith(set_description, "B.1.351") ~ "B.1.351",
    endsWith(set_description, "B.1.1.7_1") ~ "B.1.1.7",
    endsWith(set_description, "A_1") ~ "A_1",
    endsWith(set_description, "Unassigned") ~ "Unassigned",
    endsWith(set_description, "B.1.617.2_7") ~ "B.1.617.2",
    endsWith(set_description, "AY.48") ~ "AY.*",
    endsWith(set_description, "B.1.617.2_2") ~ "B.1.617.2",
    endsWith(set_description, "A_11") ~ "A_11",
    endsWith(set_description, "B.1.621") ~ "B.1.621",
    endsWith(set_description, "AY.33.2") ~ "AY.*",
    endsWith(set_description, "B.1.621.1") ~ "B.1.621",
    endsWith(set_description, "C.37") ~ "C.37",
    endsWith(set_description, "AY.20") ~ "AY.*",
    endsWith(set_description, "A_2") ~ "A_2",
    endsWith(set_description, "BA.1.15") ~ "BA.1.*",
    endsWith(set_description, "BA.1_4") ~ "BA.1.*"
  ))
head(df)


#############################
#Prepare data for ggalluvial
#############################

list <- df %>%
  group_by(lineage_set, lineage_genome, completeness, depth) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)

list$completeness <- as.numeric(list$completeness)
list$depth <- as.numeric(list$depth)

lin <- ggplot(list, aes(axis1 = lineage_genome, axis2= lineage_set, y = n)) +
  geom_alluvium(aes(fill = completeness), aes.bind=F, width = 1/4) +
  geom_stratum(width = 1/4, fill = "white", color = "black") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size =3) +
  scale_x_discrete(limits = c("genome", "hedgehog"),
                   expand = c(.05, .05)) +
  scale_fill_viridis(option = "F", limits = c(0, 100), name = "S gene completeness") +
  guides(fill=guide_legend(title="S gene completeness")) +
  labs(y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14)) + 
  ggtitle("Standard-S method")
lin

#Compare genome lineage with hedghehog lineage set only at lineage level
df$status_hedgehog <- df$lineage_genome == df$lineage_set
df$status_hedgehog <- gsub("FALSE", "Incorrect", df$status_hedgehog)
df$status_hedgehog <- gsub("TRUE", "Correct", df$status_hedgehog)
head(df)
dim(df)
df[which(df$set_description =="A_11"),13] <- "Correct" #checked hedgehog github to confirm

#hedgehog stats
pcth <- as.data.frame(table(df$status_hedgehog))
pcth$pct <- pcth$Freq/sum(pcth$Freq)*100
colnames(pcth) <- c("status", "count", "pct")
pcth$method <- c("hedgehog")

dir.create("hedgehog_Std-S")
dir.create("hedgehog_Std-S/tables")
dir.create("hedgehog_Std-S/tables/experimental_samples")
write.table(pcth, "hedgehog_Std-S/tables/experimental_samples/hedgehog_pct.tsv", sep = "\t", row.names = F, quote = F)

incorrecth <- df[which(df$status_hedgehog == "Incorrect"),]
table(incorrecth$genome_lineage)
mean(incorrecth$depth)
sd(incorrecth$depth)
mean(incorrecth$completeness)
sd(incorrecth$completeness)

incorrecth <- incorrecth[-which(incorrecth$genome_lineage == "B.1.621"),]
unique(incorrecth$genome_lineage)
table(incorrecth$genome_lineage)

dim(incorrecth)
mean(incorrecth$depth)
sd(incorrecth$depth)
mean(incorrecth$completeness)
sd(incorrecth$completeness)
max(incorrecth$completeness)
min(incorrecth$completeness)

correcth <- df[which(df$status_hedgehog == "Correct"),]
mean(correcth$depth)
sd(correcth$depth)
mean(correcth$completeness)
sd(correcth$completeness)

#Completeness boxplot
difc <- ggplot(df, aes(status_hedgehog, completeness)) + 
  geom_boxplot(alpha = 0) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14, angle = 60, hjust = 1, vjust = 0.9),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) +
  xlab("") + ylab("Completeness (%)") +
  scale_fill_viridis(option = "F", limits = c(0, 100)) +
  guides(fill=guide_legend(title="")) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  geom_jitter(size = 5, pch=21, aes(fill = completeness), alpha = 0.5) +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T, test = "wilcox.test")
difc

#Depth boxplot
difd <- ggplot(df, aes(status_hedgehog, depth)) + 
  geom_boxplot(alpha = 0) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14, angle = 60, hjust = 1, vjust = 0.9),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) +
  xlab("") + ylab("Depth (X)") +
  scale_fill_viridis(option = "F", limit =c(0, 600)) +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  guides(fill=guide_legend(title="")) +
  geom_jitter(size = 5, pch=21, aes(fill = depth), alpha = 0.5) +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T, test = "wilcox.test")
difd

write.table(df, "hedgehog_Std-S/tables/experimental_samples/std_samples_hedgehog_compare.tsv", sep = "\t", row.names = F, quote = F)

library(ggpubr)

figa <- ggarrange(lin, legend ="bottom")
figb <- ggarrange(difd, difc, ncol = 2, labels = c("B", ""))

fig1 <- ggarrange(figa, figb, ncol = 2, labels = c("A"), widths = c(0.6,0.4), heights = c(0.6,0.4), common.legend = T, legend = "bottom")
fig1

###############################
#Standard sampled read analysis
###############################

merged <- read.table("./inputs/std-S_sampling_complete_table.tsv", sep = "\t", header = T)
head(merged)
dim(merged)
unique(merged$seq)

#Get omicron samples obtained with the standard amplicon sequencing protocol
merged2 <- merged[-which(merged$seq == "Std-S"),] #remove data from seq191021

#Omicron samples

targeto <- c("CUY42-007288",
             "CUY42-007291",
             "CUY42-007286",
             "CUY42-007294",
             "CUY42-007289")

omicron <- filter(merged2, id.y %in% targeto)
unique(omicron$barcode)

lin_omicron <- omicron
unique(lin_omicron$hedgehog_precision)

lin_omicron <- lin_omicron %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "A_1") ~ "Incorrect",
    endsWith(hedgehog_precision, "BA.1.15") ~ "Correct",
    endsWith(hedgehog_precision, "A_2") ~ "Incorrect",
    endsWith(hedgehog_precision, "BA.1.17.2") ~ "Correct",
    endsWith(hedgehog_precision, "BA.2.12") ~ "Incorrect",
    endsWith(hedgehog_precision, "AY.103") ~ "Incorrect",
    endsWith(hedgehog_precision, "A.23.1") ~ "Incorrect",
    endsWith(hedgehog_precision, "BA.1_4") ~ "Correct",
    endsWith(hedgehog_precision, "BA.1.1_3") ~ "Correct"
  ))
head(lin_omicron)

lin_omicron$sampling <- as.numeric(as.character(lin_omicron$sampling))
lin_omicron$depth <- as.numeric(as.character(lin_omicron$depth))
lin_omicron$completeness <- as.numeric(as.character(lin_omicron$completeness))


omicron1 <- ggplot(lin_omicron, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10), 
                     legend.title = element_text(size = 18), 
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="hedgehog")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Omicron")
omicron1 

omicron2 <-  ggplot(lin_omicron, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.1, size = 3) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10), 
                     legend.title = element_text(size = 18), 
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="hedgehog")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
omicron2


#remove barcodes from seq3
merged3 <- merged[-which(merged$seq == "omicron"),]


#alpha

targeta <- c("CUY17-003849",
             "CUY17-003866",
             "CUY17-003867")

alpha <- filter(merged3, id.y %in% targeta)
unique(alpha$barcode)

#beta

targetb <- c("CUY16-003804",
             "CUY16-003803",
             "CUY16-003800",
             "CUY16-003802",
             "CUY16-003801")
beta <- filter(merged3, id.y %in% targetb)
unique(beta$barcode)


#Gamma
targetg <- c("CUY17-003939",
             "CUY12-002713",
             "CUY12-002716",
             "CUY12-002721",
             "CUY12-002725")

gamma <- filter(merged3, id.y %in% targetg)
unique(gamma$barcode)


#Delta

targetd <- c("CUY26-005235",
             "CUY17-003909",
             "CUY24-004954",
             "CUY23-004816",
             "CUY16-003792",
             "CUY24-004876",
             "CUY17-003893",
             "CUY17-003901",
             "CUY17-003892",
             "CUY17-003910")

delta <- filter(merged, id.y %in% targetd)
dim(delta)
unique(delta$barcode)


#P.6
targetp <- c("CUY5-000547",
             "CUY16-003498",
             "CUY16-003501",
             "CUY16-003536",
             "CUY16-003538")
p6 <- filter(merged3, id.y %in% targetp)   
unique(p6$barcode)


#Mu

targetm <- c("CUY21-004441",
             "CUY18-004026",
             "CUY18-004019",
             "CUY18-004022",
             "CUY18-004023")

mu <- filter(merged, id.y %in% targetm)   
unique(mu$barcode)

#Lambda

targetl <- c("CUY19-004183",
             "CUY6-001233",
             "CUY8-001646")
lambda <- filter(merged3, id.y %in% targetl)   
unique(lambda$barcode)

#P.2
targetpp <- c("CUY4-000486",
              "CUY16-003496",
              "CUY1-000006")
p2 <- filter(merged3, id.y %in% targetpp)   
unique(p2$barcode)



##############
###ALPHA#######
##############

lin_alfa <- alpha
unique(lin_alfa$hedgehog_precision)
unique(lin_alfa$barcode)

lin_alfa$PANGO_lineage <- ifelse(grepl("B.1.1.7_1", lin_alfa$hedgehog_precision, ignore.case = T), "Correct",
                                 ifelse(grepl("", lin_alfa$hedgehog_precision, ignore.case = T), "Incorrect", "Incorrect"))

lin_alfa$sampling <- as.numeric(as.character(lin_alfa$sampling))
lin_alfa$depth <- as.numeric(as.character(lin_alfa$depth))
lin_alfa$completeness <- as.numeric(as.character(lin_alfa$completeness))

a <- ggplot(lin_alfa, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Alpha") + labs(color='PANGO lineage')
a 

aa <-  ggplot(lin_alfa, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") +
  labs(color='PANGO lineage')
aa


##############
###BETA#######
##############

lin_beta <- beta
unique(lin_beta$barcode)
unique(lin_beta$hedgehog_precision)

lin_beta <- lin_beta %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "B.1.351") ~ "Correct",
    endsWith(hedgehog_precision, "A_23") ~ "Incorrect", 
    endsWith(hedgehog_precision, "A_1") ~ "Incorrect"))
head(lin_beta)


lin_beta$sampling <- as.numeric(as.character(lin_beta$sampling))
lin_beta$depth <- as.numeric(as.character(lin_beta$depth))
lin_beta$completeness <- as.numeric(as.character(lin_beta$completeness))

b <- ggplot(lin_beta, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Beta") + labs(color='PANGO lineage')
b 

bb <-  ggplot(lin_beta, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
bb

##############
###GAMMA######
##############


lin_gamma <- gamma
unique(lin_gamma$barcode)
unique(lin_gamma$hedgehog_precision)


lin_gamma <- lin_gamma %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "P.1_1") ~ "Correct"))
head(lin_gamma)


lin_gamma$sampling <- as.numeric(as.character(lin_gamma$sampling))
lin_gamma$depth <- as.numeric(as.character(lin_gamma$depth))
lin_gamma$completeness <- as.numeric(as.character(lin_gamma$completeness))


c <- ggplot(lin_gamma, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Gamma") + labs(color='PANGO lineage')
c

cc <- ggplot(lin_gamma, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
cc

##############
###DELTA######
##############

lin_delta <- delta
unique(lin_delta$barcode)
unique(lin_delta$hedgehog_precision)

lin_delta[is.na(lin_delta)] <- "Unassigned"

lin_delta <- lin_delta %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "B.1.617.2_2") ~ "Correct",
    endsWith(hedgehog_precision, "A_11") ~ "Correct",
    endsWith(hedgehog_precision, "AY.126") ~ "Correct",
    endsWith(hedgehog_precision, "AY.103") ~ "Correct",
    endsWith(hedgehog_precision, "Unassigned") ~ "Incorrect",
    endsWith(hedgehog_precision, "A_1") ~ "Incorrect",
    endsWith(hedgehog_precision, "A_2") ~ "Incorrect",
    endsWith(hedgehog_precision, "B.1.551") ~ "Incorrect",
    endsWith(hedgehog_precision, "B.1.1.523") ~ "Incorrect",
    endsWith(hedgehog_precision, "AY.20") ~ "Correct",
    endsWith(hedgehog_precision, "B.1.617.2_1") ~ "Correct",
    endsWith(hedgehog_precision, "AY.30") ~ "Correct",
    endsWith(hedgehog_precision, "B.1.617.2_7") ~ "Correct"
  ))
head(lin_delta)



lin_delta$sampling <- as.numeric(as.character(lin_delta$sampling))
lin_delta$depth <- as.numeric(as.character(lin_delta$depth))
lin_delta$completeness <- as.numeric(as.character(lin_delta$completeness))

d <- ggplot(lin_delta, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Delta") + labs(color='PANGO lineage')
d 

dd <- ggplot(lin_delta, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
dd


#############
###P.6#######
#############

lin_p <- p6
unique(lin_p$hedgehog_precision)
unique(lin_p$barcode)

lin_p <- lin_p %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "P.6") ~ "Correct",
    endsWith(hedgehog_precision, "A_1") ~ "Incorrect"))
head(lin_p)

lin_p$sampling <- as.numeric(as.character(lin_p$sampling))
lin_p$depth <- as.numeric(as.character(lin_p$depth))
lin_p$completeness <- as.numeric(as.character(lin_p$completeness))

e <- ggplot(lin_p, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("P.6") + labs(color='PANGO lineage')
e

ee <- ggplot(lin_p, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
ee


#############
###MU########
#############

lin_mu <- mu
unique(lin_mu$barcode)
unique(lin_mu$hedgehog_precision)

lin_mu <- lin_mu %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "B.1.621") ~ "Correct",
    endsWith(hedgehog_precision, "B.1.617.2_2") ~ "Incorrect",
    endsWith(hedgehog_precision, "AY.33.2") ~ "Incorrect",
    endsWith(hedgehog_precision, "A.30") ~ "Incorrect",
    endsWith(hedgehog_precision, "B.1.621.1") ~ "Correct"
  ))
head(lin_mu)

lin_mu$sampling <- as.numeric(as.character(lin_mu$sampling))
lin_mu$depth <- as.numeric(as.character(lin_mu$depth))
lin_mu$completeness <- as.numeric(as.character(lin_mu$completeness))



f <- ggplot(lin_mu, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Mu") + labs(color='PANGO lineage')
f

ff <- ggplot(lin_mu, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
ff


#############
###LAMBDA####
#############

lin_lambda <- lambda
unique(lin_lambda$barcode)
unique(lin_lambda$hedgehog_precision)

lin_lambda <- lin_lambda %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog, "C.37") ~ "Correct",
    endsWith(hedgehog, "C.37.1") ~ "Correct"))
head(lin_lambda)

lin_lambda$sampling <- as.numeric(as.character(lin_lambda$sampling))
lin_lambda$depth <- as.numeric(as.character(lin_lambda$depth))
lin_lambda$completeness <- as.numeric(as.character(lin_lambda$completeness))


h <- ggplot(lin_lambda, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("Lambda") + labs(color='PANGO lineage')
h

hh <- ggplot(lin_lambda, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
hh

#############
###P.2#######
#############

lin_p2 <- p2
unique(p2$barcode)
unique(lin_p2$hedgehog_precision)

lin_p2 <- lin_p2 %>%
  mutate(PANGO_lineage = case_when(
    endsWith(hedgehog_precision, "P.2") ~ "Correct",
    endsWith(hedgehog_precision, "A_23") ~ "Incorrect"))
head(lin_p2)


lin_p2$sampling <- as.numeric(as.character(lin_p2$sampling))
lin_p2$depth <- as.numeric(as.character(lin_p2$depth))
lin_p2$completeness <- as.numeric(as.character(lin_p2$completeness))

i <- ggplot(lin_p2, aes(x=sampling, y=depth, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 500)) +  guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("P.2") + labs(color='PANGO lineage')
i


ii <- ggplot(lin_p2, aes(x=sampling, y=completeness, color = PANGO_lineage)) + geom_point(alpha = 0.5, size = 3) +
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 10),
                     legend.title = element_text(size = 18),
                     legend.text = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 100)) + guides(color=guide_legend(title="lineage set")) +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  xlab("") + ylab("") + ggtitle("") + labs(color='PANGO lineage')
ii

######################
##Supplementary S1
######################

library(ggpubr)
library(grid)
library(gridExtra)

f1a <- ggarrange(a, b, c, d, omicron1,  ncol = 5, nrow = 1, common.legend = T,
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f1a <- annotate_figure(f1a, left = textGrob("Average sequencing depth (X)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("", gp = gpar(cex = 1.3)))
f1a 


f1b <- ggarrange(aa, bb, cc, dd, omicron2, ncol = 5, nrow = 1, legend = "none",
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom"))

f1b
f1b <- annotate_figure(f1b, left = textGrob("% S gene completeness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Sampled reads", gp = gpar(cex = 1.3)))
f1b


S1 <- ggarrange(f1a, f1b, ncol = 1, nrow = 2, common.legend = T,
                align = "hv",
                labels = c("A", "B"),
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
S1



dir.create("hedgehog_Std-S/supplementary1")
png('./hedgehog_Std-S/supplementary1/voc_depth_completness.png', res = 600, height = 25, width = 40, units = 'cm')
S1
dev.off()


#####################
##Supplementary S2
#####################

f1a <- ggarrange(f, h, i, e,  ncol = 4, nrow = 1, common.legend = T,
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f1a <- annotate_figure(f1a, left = textGrob("Average sequencing depth (X)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Sampled reads", gp = gpar(cex = 1.3)))
f1a 

f1b <- ggarrange(ff, hh, ii, ee, ncol = 4, nrow = 1, common.legend = T,
                 align = "hv",
                 font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom"))

f1b

f1b <- annotate_figure(f1b, left = textGrob("% S gene completeness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Sampled reads", gp = gpar(cex = 1.3)))
f1b


S2 <- ggarrange(f1a, f1b, ncol = 1, nrow = 2, common.legend = T,
                align = "hv",
                labels = c("A", "B"),
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
S2

dir.create("hedgehog_Std-S/supplementary2")
png('./hedgehog_Std-S/supplementary2/novocvoi_depth_completness.png', res = 600, height = 25, width = 40, units = 'cm')
S2
dev.off()


##################
#Sampling analysis
##################

tab <- as.data.frame(rbind(lin_alfa, lin_beta, lin_gamma, lin_delta, lin_omicron, lin_mu, lin_lambda, lin_p2, lin_p))
head(tab)
dim(tab)
tab$completeness <- as.numeric(as.character(tab$completeness))
tab$depth <- as.numeric(as.character(tab$depth))
head(tab)
dim(tab)
colnames(tab)[12] <- c("status")

dir.create("hedgehog_Std-S/tables/sampling_std_hedgehog")
write.table(tab, "hedgehog_Std-S/tables/sampling_std_hedgehog/hedgehog_std_sampling.tsv", sep = "\t", row.names = F, quote = F)


failed <- tab[which(tab$status == "Incorrect"),]
head(failed)
mean_faield_depth <- as.data.frame(mean(failed$depth))
sd_failed_depth <- as.data.frame(sd(failed$depth))
failed_stats_depth <- as.data.frame(cbind(mean_faield_depth, sd_failed_depth))
failed_stats_depth$maximum <- max(failed$depth)
failed_stats_depth$minimum <- min(failed$depth)
rownames(failed_stats_depth) <- c("Incorrect")
colnames(failed_stats_depth) <- c("mean", "SD", "max", "min")
failed_stats_depth <- round(failed_stats_depth, 2)

mean_failed_completness <- as.data.frame(mean(failed$completeness))
sd_failed_completeness <- as.data.frame(sd(failed$completeness))
failed_stats_completeness <- as.data.frame(cbind(mean_failed_completness, sd_failed_completeness))
failed_stats_completeness$maximum <- max(failed$completeness)
failed_stats_completeness$minimum <- min(failed$completeness)
rownames(failed_stats_completeness) <- c("Incorrect")
colnames(failed_stats_completeness) <- c("mean", "SD", "max", "min")
failed_stats_completeness <- round(failed_stats_completeness, 2)

pass <- tab[which(tab$status == "Correct"),]
mean_pass_depth <- as.data.frame(mean(pass$depth))
sd_pass_depth <- as.data.frame(sd(pass$depth))
pass_stats_depth <- as.data.frame(cbind(mean_pass_depth, sd_pass_depth))
pass_stats_depth$maximum <- max(pass$depth)
pass_stats_depth$minimum <- min(pass$depth)
rownames(pass_stats_depth) <- c("Correct")
colnames(pass_stats_depth) <- c("mean", "SD", "max", "min")
pass_stats_depth <- round(pass_stats_depth, 2)

mean_pass_completness <- as.data.frame(mean(pass$completeness))
sd_pass_completeness <- as.data.frame(sd(pass$completeness))
pass_stats_completeness <- as.data.frame(cbind(mean_pass_completness, sd_pass_completeness))
pass_stats_completeness$maximum <- max(pass$completeness)
pass_stats_completeness$minimum <- min(pass$completeness)
rownames(pass_stats_completeness) <- c("Correct")
colnames(pass_stats_completeness) <- c("mean", "SD", "max", "min")
pass_stats_completeness <- round(pass_stats_completeness, 2)

summary_pass_depth <- as.data.frame(rbind(pass_stats_depth, failed_stats_depth))
summary_pass_depth$status <- row.names(summary_pass_depth)
summary_completeness <- as.data.frame(rbind(pass_stats_completeness, failed_stats_completeness))
summary_completeness$status <- row.names(summary_completeness)

write.table(summary_pass_depth, "hedgehog_Std-S/tables/sampling_std_hedgehog/hedgehog_depth_sampling.tsv", sep = "\t", row.names = F, quote = F)
write.table(summary_completeness, "hedgehog_Std-S/tables/sampling_std_hedgehog/hedgehog_completeness_sampling.tsv", sep = "\t", row.names = F, quote = F)

#Boxplot depth sampling
depth_g <- ggplot(tab, aes(status, depth, fill = status)) + 
  geom_boxplot(alpha = 0.5) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.position = "none") +
  guides(fill=guide_legend(title="Status")) +
  xlab("") + ylab("Depth (X)") + 
  scale_fill_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  #scale_fill_viridis(option = "F", discrete = T, direction = -1) +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T)
depth_g

#Boxplot completeness sampling 
comp_g <- ggplot(tab, aes(status, completeness, fill = status)) + 
  geom_boxplot(alpha = 0.5) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.position = "none") +
  guides(fill=guide_legend(title="Status")) +
  xlab("") + ylab("Completeness (%)") + 
  scale_fill_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059")) +
  #scale_fill_viridis(option = "F", discrete = T, direction = -1) +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T)
comp_g

stats <- ggarrange(depth_g, comp_g, ncol = 2, nrow = 1, common.legend = T,  legend = "none",
                   font.label = list(size = 18, color = "black", face = "bold", family = NULL, position = "top"))
stats

#Dotplot depth and completeness
cor <- ggplot(tab, aes(x=depth, y=completeness, color = status)) + geom_point(alpha = 0.3, size = 5) +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     legend.position = "right",
                     panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(size = 14, face = "bold"),
                     axis.title = element_text(size = 14, face = "bold"),
                     axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  geom_vline(xintercept =312, linetype = "dashed")+
  geom_hline(yintercept = 90, linetype = "dashed") +
  annotate("text", x=390, y=91, label= "min completeness") +
  annotate("text", x=313, y=97, label= "min depth", angle = 90) +
  xlab("Average depth (X)") + ylab("Completeness (%)") + 
  labs(color='lineage set status') +
  scale_color_manual(values = c("Correct" = "#fcd5ce", "Incorrect" = "#9e0059"))
cor


fig2 <- ggarrange(stats, cor, ncol = 2, nrow = 1, common.legend = F,
                  labels = c("C", "D"), font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))
fig2


FIG1f <- ggarrange(fig1, fig2, nrow = 2, widths = c(0.5, 0.5), heights = c(0.5, 0.5))
FIG1f 

dir.create("hedgehog_Std-S/figure1")
png('./hedgehog_Std-S/figure1/figure1.png', res = 600, height = 35, width = 40, units = 'cm')
FIG1f 
dev.off()


d <- tab %>% filter(completeness > 90)
dd <- tab %>% filter(depth > 312)
table(dd$status)
mean(dd$sampling)
min(dd$sampling)
max(dd$sampling)

t <- as.data.frame(table(dd$status))
t$pct <- t$Freq/sum(t$Freq)*100
colnames(t) <- c("status", "count", "pct")
t
write.table(t, "hedgehog_Std-S/tables/sampling_std_hedgehog/quadrant_depth_compl.tsv", sep = "\t", row.names = F, quote = F)
