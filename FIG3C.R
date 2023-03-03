library(ggplot2)
library(stringr)
library(ggsignif)
library(dplyr)
library(ggalluvial)
library(tidyr)

df <- read.table("./inputs/gisaid_omicron_table.tsv", sep = "\t", header = T)
head(df)
colnames(df)[2] <- c("genome_lineage")

dir.create("gisaid_omicron/")

sp1 <- as.data.frame(str_split_fixed(df$genome_lineage, "[.]",  3))
a <- as.data.frame(str_c(sp1$V1,".", sp1$V2))
sp2 <- as.data.frame(str_split_fixed(df$hedgehog_set_hash, "[.]",  3))
b <- as.data.frame(str_c(sp2$V1,".", sp2$V2))

new <- as.data.frame(cbind(a, b))
colnames(new) <- c("genome", "hedgehog")
head(new)

df <- as.data.frame(cbind(df, new))
head(df)
dim(df)

df$status_hedgehog <- df$genome == df$hedgehog
table(df$status_hedgehog)

df$status_hedgehog <- gsub("TRUE", "Correct", df$status_hedgehog)
df$status_hedgehog <- gsub("FALSE", "Incorrect", df$status_hedgehog)

head(df)

hg <- read.csv("./inputs/set_names.95.csv", header = T)
head(hg)

hg[which(hg$set_description == "A_1"),]

target <- unique(df$hedgehog_precision)
hdg <- hg
head(hdg)
hdg <- filter(hdg, set_description %in% target)
colnames(hdg)[5] <- c("hedgehog_precision")
hdg <- select(hdg, hedgehog_precision, lineages)

df <- merge(df, hdg, by = "hedgehog_precision")
head(df)
dim(df)

#####################
###Lineage level only
#####################

a <- df[which(df$hedgehog_set_hash == "A"),]
unique(a$hedgehog_precision)

df[which(df$hedgehog_precision =="A_14"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_17"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_2"),9] <- "Incorrect"
df[which(df$hedgehog_precision =="A_1"),9] <- "Incorrect"
df[which(df$hedgehog_precision =="A_4"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_3"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_12"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_15"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_25"),9] <- "Correct"
df[which(df$hedgehog_precision =="A_9"),9] <- "Correct"

dim(df)
head(df)

pattern <- df$genome_lineage
dff <- df %>% 
  mutate(across(lineages, ~case_when(str_detect(., pattern) ~ str_extract(., pattern)), .names = "new_col{col}")) %>% 
  unite(lineage_set, starts_with('new'), na.rm = TRUE, sep = '\t')
head(dff)
dff[dff == ''] <- "unmatched"

dff$status_lineageset <- ifelse(grepl("unmatched", dff$lineage_set, ignore.case = T), "Incorrect",
                                ifelse(grepl("", dff$lineage_set, ignore.case = T), "Correct", "Correct"))

unique(df$hedgehog)
dff <- dff %>%
  mutate(hedgehog2 = case_when(
    endsWith(hedgehog, "BA.2") ~ "BA.2.*",
    endsWith(hedgehog, "BA.1") ~ "BA.1.*",
    endsWith(hedgehog, "A.") ~ "A",
    endsWith(hedgehog, "BD.1") ~ "BD.1.*",
    endsWith(hedgehog, "B.1") ~ "B.1.*",
    endsWith(hedgehog, "AY.4") ~ "AY.*", 
    endsWith(hedgehog, "AY.43") ~ "AY.*", 
    endsWith(hedgehog, "XAL.") ~ "XAL",
    endsWith(hedgehog, "XAR.") ~ "XAR",
    endsWith(hedgehog, "BC.1") ~ "BC.1.*",
    endsWith(hedgehog, "AY.53") ~ "AY.*",
    endsWith(hedgehog, "XAM.") ~ "XAM.*",
    endsWith(hedgehog, "BA.5") ~ "BA.5.*",
    endsWith(hedgehog, "BA.4") ~ "BA.4.*",
    endsWith(hedgehog, "P.1") ~ "P.1.*",
    endsWith(hedgehog, "DE.1") ~ "DE.1.*",
    endsWith(hedgehog, "BF.6") ~ "BF.6.*"
  ))

unique(df$genome)
dff <- dff %>%
  mutate(genome2 = case_when(
    endsWith(genome, "BA.2") ~ "BA.2.*",
    endsWith(genome, "BA.1") ~ "BA.1.*",
    endsWith(genome, "XG.") ~ "XG.*",
    endsWith(genome, "BC.1") ~ "BC.1.*",
    endsWith(genome, "BD.1") ~ "BD.1.*",
    endsWith(genome, "XS.") ~ "XS.*", 
    endsWith(genome, "XE.") ~ "XE.*", 
    endsWith(genome, "XR.") ~ "XR.*",
    endsWith(genome, "XQ.") ~ "XQ.*",
    endsWith(genome, "XM.") ~ "XM.*",
    endsWith(genome, "BA.5") ~ "BA.5.*",
    endsWith(genome, "BA.4") ~ "BA.4.*",
    endsWith(genome, "XAA.") ~ "XAA.*",
    endsWith(genome, "BE.1") ~ "BE.1.*",
    endsWith(genome, "BE.3") ~ "BE.3.*",
    endsWith(genome, "Unassigned.") ~ "Unassigned",
    endsWith(genome, "BF.5") ~ "BF.5.*",
    endsWith(genome, "BC.2") ~ "BC.2.*",
    endsWith(genome, "BA.3") ~ "BA.3.*"
  ))


dff <- dff[-which(dff$genome_lineage == "Unassigned"),]

table(dff$status_hedgehog)
table(dff$status_lineageset)

#correct vs. incorrect reduced nomeclature
th <- select(dff, genome_lineage, hedgehog2, status_hedgehog, pct) 
th$method <- c("hedgehog")
colnames(th)[2] <- c("lineage")
colnames(th)[3] <- c("status")
head(th)


e <- ggplot(th, aes(status, pct)) + 
  geom_boxplot(alpha = 0.5) +
  theme(text = element_text(size = 25),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        title = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        strip.background =element_rect(fill="white"),
        axis.title = element_text(size = 20, face = "bold")) +
  xlab("") + ylab("S gene completeness (%)") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 120, by = 25)) +
  guides(fill=guide_legend(title="")) +
  ggtitle("", subtitle = "") +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T, test = "wilcox.test", y_position = c(105),textsize=6)
e

#correct vs. incorrect complete matched
head(dff)
thh <- select(dff, genome_lineage, hedgehog_precision, status_lineageset, pct) 
thh$method <- c("hedgehog")
colnames(thh)[2] <- c("exact lineage")
colnames(thh)[3] <- c("status")
head(thh)

f <- ggplot(th, aes(status, pct)) + 
  geom_boxplot(alpha = 0.5) +
  theme(text = element_text(size = 25),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        title = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        strip.background =element_rect(fill="white"),
        axis.title = element_text(size = 12)) +
  xlab("") + ylab("S gene completeness (%)") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 105, by = 25)) +
  guides(fill=guide_legend(title="")) +
  ggtitle("", subtitle = "") +
  geom_signif(comparisons = list(c("Correct", "Incorrect")), map_signif_level = T, test = "wilcox.test", y_position = c(105),textsize=6)
f

png("./gisaid_omicron/figs/hedgehog_precise.png", res = 600, height = 25, width = 40, units = "cm")
f
dev.off()

###########################################################################################################################################
#Main Omicron lineages 

target <- c("BA.1", "BA.2", "BA.4", "BA.5")
dt <- filter(dff, genome %in% target)
dim(dt)
unique(dt$genome)

tab <- dt
head(tab)
dim(tab)

hedge <- select(tab, accession, genome_lineage, genome2, hedgehog_set_hash, hedgehog, status_hedgehog)
head(hedge)

hedge <- hedge %>%
  group_by(genome_lineage, genome2, hedgehog_set_hash, status_hedgehog) %>%
  count(status = status_hedgehog)
hedge <- as.data.frame(hedge)

#Main figure 3, B
p1 <- ggplot(hedge, aes(x = genome2, y = n, fill =status_hedgehog))
s4 <- p1 +  geom_bar(alpha = 0.8, stat="identity", position="fill") +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.8),
        legend.text=element_text(size=12),
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(nrow=1, ncol = 2)) + 
  scale_fill_manual(values = c("Correct"="#f28482", "Incorrect"="#f9dcc4")) +
  xlab("") + ylab("") + labs(fill = "Status")
s4


#Supplementary figure 
hedge <- select(tab, accession, genome_lineage, genome2, status_lineageset, hedgehog_set_hash, pct)
head(hedge)
pct <- as.data.frame(table(hedge$genome2, hedge$status_lineageset))
pct

inc <- hedge[which(hedge$status_lineageset == "Incorrect"),]
mean(inc$pct)
sd(inc$pct)

cor <- hedge[which(hedge$status_lineageset == "Correct"),]
mean(cor$pct)
sd(cor$pct)

hedge <- hedge %>%
  group_by(genome_lineage, genome2, status_lineageset) %>%
  count(status = status_lineageset)
hedge <- as.data.frame(hedge)



p1 <- ggplot(hedge, aes(x = genome_lineage, y = n, fill = status_lineageset))
s5 <- p1 +  geom_bar(alpha = 0.8, stat="identity", position="fill") +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_text(size = 14),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.8),
        legend.text=element_text(size=12),
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(nrow=1, ncol = 2)) + 
  facet_grid(~genome2, scales = "free_x") +
  ggtitle("S gene PANGO lineage assigment") +
  scale_fill_manual(values = c("Correct"="#f28482", "Incorrect"="#f9dcc4")) +
  xlab("Genome PANGO lineage") + ylab("S gene based PANGO lineage assignment") + labs(fill = "Status")
s5


library(ggpubr)


fig <- ggarrange(s5, nrow = 2, labels = c("A", "B"), widths = c(0.6, 0.4), heights = c(0.6, 0.4))
fig

dir.create("gisaid_omicron")
dir.create("gisaid_omicron/figs")
png("./gisaid_omicron/figs/supp_hedgehog_exact_lineage.png", res = 600, height = 30, width = 40, units = "cm")
fig
dev.off()

####################
#Main figure 3
####################

ba1 <- dff[grepl("BA.1.*", dff$genome2),]
a <- as.data.frame(table(ba1$status_hedgehog))
a$pct <- a$Freq/sum(a$Freq)*100
a$method <- c("hedgehog")
table(ba1$status_usher)


ba2 <- dff[grepl("BA.2.*", dff$genome2),]
a <- as.data.frame(table(ba2$status_hedgehog))
a$pct <- a$Freq/sum(a$Freq)*100
a$method <- c("hedgehog")


ba4 <- dff[grepl("BA.4.*", dff$genome2),]
table(ba4$status_hedgehog)
table(ba4$status_usher)
a <- as.data.frame(table(ba4$status_hedgehog))
a$pct <- a$Freq/sum(a$Freq)*100
a$method <- c("hedgehog")

ba5 <- dff[grepl("BA.5.*", dff$genome2),]
table(ba5$status_hedgehog)
a <- as.data.frame(table(ba5$status_hedgehog))
a$pct <- a$Freq/sum(a$Freq)*100
a$method <- c("hedgehog")

head(dff)
s <- dff
s <- select(s, accession, genome_lineage,  hedgehog_set_hash, hedgehog_precision, 
            genome2, hedgehog2,  status_hedgehog, status_lineageset, lineages)
head(s)
colnames(s) <- c("GISAID_accession", "GISAID_lineage",  "hedgehog_precision", "hedgehog_set-description", 
                 "compre_genome", "compare_hedgehog", "status_hedgehog", "status_lineageset", "lineages")

dir.create("gisaid_omicron")
dir.create("gisaid_omicron/tables")
write.table(s, "./gisaid_omicron/tables/supplemetary_gisaid_omicron.tsv", sep = "\t", row.names = F, quote = F)

#################################################################################################################
target <- c("BA.1", "BA.2", "BA.4", "BA.5")
dt <- filter(dff, genome %in% target)
dim(dt)
unique(dt$genome)

set.seed(50)
dff <- dt[sample(nrow(dt), 300), ]
dim(dff)
head(dff)
unique(dff$hedgehog_precision)

list <- dff %>%
  group_by(genome, genome_lineage, genome2, hedgehog_precision, status_hedgehog, status_lineageset) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)
unique(list$genome)

#Main figure 3 A
lin3 <- ggplot(list, aes(axis1 = genome2, axis2= hedgehog_precision)) +
  geom_alluvium(aes(fill = status_hedgehog), aes.bind=F, width = 1/4, alpha = 0.8) +
  geom_stratum(width = 1/4, fill = "white", color = "black") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size =4) +
  scale_x_discrete(limits = c("GISAID", "hedgehog"),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = c("Correct"="#f28482", "Incorrect"="#f9dcc4")) +
  labs(y = "") + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12), 
        plot.title = element_text(size = 12)) + guides(fill=guide_legend(title="")) +
  ggtitle("", subtitle = "")
lin3

#####################################################################################################################

tab <- dt
head(tab)
dim(tab)

hedge <- select(tab, accession, genome_lineage, genome2, hedgehog_set_hash, hedgehog, status_hedgehog)
head(hedge)
hedge_sum <- table(hedge$genome2, hedge$status_hedgehog)
hedge_sum <- as.data.frame(hedge_sum)
head(hedge_sum)
colnames(hedge_sum) <- c("Genome", "Status", "Count")

p1 <- ggplot(hedge_sum, aes(x = Genome, y = Count, fill = Status))
s1 <- p1 +  geom_bar(alpha = 0.8, stat="identity", position="fill") +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text=element_text(size=12),
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(nrow=1, ncol = 2)) + ggtitle("") +
  scale_fill_manual(values = c("Correct"="#f28482", "Incorrect"="#f9dcc4")) +
  xlab("") + ylab("Correct vs. Incorrect assignments") + labs(fill = "Status")
s1

fig <- ggarrange(lin3, s1, ncol = 2, labels = c("A", "B"), common.legend = T, legend = "bottom")
fig

dir.create("gisaid_omicron/figure3")
png("./gisaid_omicron/figure3/figure3.png", res = 600, height = 25, width = 40, units = "cm")
fig
dev.off()
