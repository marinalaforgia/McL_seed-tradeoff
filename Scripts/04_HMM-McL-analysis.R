#### HMMs McLaughlin Manuscript ####

rm(list = ls())

calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

#### Load libraries ####
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(rdacca.hp)
library(Metrics)
library(ggcorrplot)
library(taxize)
library(WorldFlora)
library(rotl)
library(U.PhyloMaker)
library(phytools)
library(caper)
library(ape)

#devtools::install_github("jinyizju/V. PhyloMaker2")


#### Read in data ####
mcl.df.q <- readRDS("Intermediates/McL-species-boot.RDS") # HMM estimates
metadata <- read.csv("Data/HMM-meta-mcl.csv")
# seed.species <- read.csv("/Users/marina.laforgia/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/20211001_Full-Species-List.csv")
seed.traits <- read.csv("/Users/marina.laforgia/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/20230801_Seed-Traits_clean_site.csv")

mcl.df.q <- mcl.df.q %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

# update names for merger with trait data
mcl.df.q$Species_Name <- recode_factor(mcl.df.q$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

mcl.df.q <- merge(mcl.df.q, metadata, by = "Species_Name")

mcl.df.q[!mcl.df.q$Species_Name %in% metadata$Species_Name,] # none

#### Filter data ####
mcl.df.q <- filter(mcl.df.q, n.plots >= 50) # 155 species to 62 species

mcl.df.q <- mcl.df.q[mcl.df.q$s.boot.sd < 0.1 & mcl.df.q$c.boot.sd < 0.1,] # 2 species lost

mcl.df.q <- mcl.df.q[mcl.df.q$perc.con > 0.5,] # zero species lost

mcl.df.q$FunGroup <- recode_factor(mcl.df.q$FunGroup, 'Exotic Forb' = "Non-native forb", 'Exotic Grass' = "Non-native grass", 'Native Forb' = "Native forb", 'Native Grass' = "Native grass")

colnames(mcl.df.q)[2:6] <- c("p0", "g", "c", "s", "r")

#### Seed trait prep ####

# some species were collected from two sites for larger cross-site study, use those from McLaughlin
seed.traits2 <- filter(seed.traits, grepl("MCL", seed.traits$site, fixed = TRUE))
seed.traits2 <- filter(seed.traits2, !grepl("CP", site))
seed.traits <- filter(seed.traits, Species %in% mcl.df.q$Species_Name, !grepl("CP", site))
seed.traits <- unique(rbind(seed.traits, seed.traits2))
rm(seed.traits2)

# create functional group factor
seed.traits$nat.inv <- recode_factor(seed.traits$nat.inv, invasive = "Non-native", native = "Native")
seed.traits$FunGroup <- paste(seed.traits$nat.inv, seed.traits$group, sep = " ")


# new mass is persistent appendages only (indehiscent fruits), chem mass includes achenes, mericarps, and nutlets where removing fruit would be nearly impossible
pers <- c("MICCAL", "MICDOU", "EROCIC", "ATHPUS", "THYCUR")
pers <- c(pers, seed.traits[seed.traits$group == "grass",]$new.code)
seed.traits$new.mass <- ifelse(seed.traits$new.code %in% pers, seed.traits$morph.mass.mg, seed.traits$chem.mass.mg) 

# Normalize seed.traits
#hist(log(seed.traits$wing.loading))
seed.traits$wing.loading <- log(seed.traits$wing.loading)

#hist(log(seed.traits$coat.perm.perc))
seed.traits$coat.perm.perc <- log(seed.traits$coat.perm.perc)

#hist(log(seed.traits$morph.mass.mg))
seed.traits$morph.mass.mg <- log(seed.traits$morph.mass.mg)

hist(log(seed.traits$chem.mass.mg))
seed.traits$chem.mass.mg <- log(seed.traits$chem.mass.mg)

# hist(log(seed.traits$new.mass))
seed.traits$new.mass <- log(seed.traits$new.mass)

#hist(log(seed.traits$size.mm))
seed.traits$size.mm <- log(seed.traits$size.mm)

mcl.df.q.all <- merge(mcl.df.q, seed.traits, by.x = c("Species_Name", "FunGroup", "new.code"), by.y = c("Species", "FunGroup", "new.code"), all.x = T) # for tree

mcl.df.q.trait <- merge(mcl.df.q, seed.traits, by.x = c("Species_Name", "FunGroup", "new.code"), by.y = c("Species", "FunGroup", "new.code"), all = F)

missing <- mcl.df.q[!mcl.df.q$Species_Name %in% mcl.df.q.trait$Species_Name,]

# 10 Species we lack traits for: 
## Native forbs: Acmispon brachycarpus, Castilleja rubicundula, Croton setiger, Eriogonum vimineum, Galium aparine, Leptosiphon bicolor, Navarretia jepsonii, Triphysaria eriantha
## Non-native-forbs: Lactuca serriola,Cerastium glomeratum

#### Colors for figs ####
fungroup_cols <- c(
  "Non-native forb" = "#d95f02",
  "Native forb" = "#7570b3",
  "Non-native grass" = "#1b9e77",
  "Native grass" = "#146c54"
)
 
fungroup_shapes <- c(
  "Non-native forb" = 17,
  "Native forb" = 19,
  "Non-native grass" = 15,
  "Native grass" = 23
)
  
#### Figure 1: Trade-off ####
##### Figure ####

# prep figure to graph CESO datapoint on top (makes it easier to see with the seed pic)
fig1.df <- mcl.df.q
fig1.df$s1 <- NA
fig1.df$c1 <- NA
fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$s1 <- fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$s
fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$c1 <- fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$c

fig1.a <- ggplot(fig1.df, aes(x = s, y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(inherit.aes = F, aes(x = s, y = c), method = "lm", col = "black", se = F) + 
  geom_point(aes(col = FunGroup, shape = FunGroup, fill = FunGroup), size = 3) +
  geom_point(aes(x = s1, y = c1, col = FunGroup, shape = FunGroup, fill = FunGroup), size = 3) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.position = c(.75,.88),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
  ) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "Probability of temporal dispersal",
       y = "Probability of spatial dispersal")

mcl.df.q$FunGroup <- factor(mcl.df.q$FunGroup, 
                            levels = c("Non-native grass", "Non-native forb", "Native forb", "Native grass"))

letters_df <- data.frame(
  FunGroup = c("Non-native grass", "Non-native forb", "Native forb"),
  label = c("a", "ab", "b")
)

fig1.b <- ggplot(mcl.df.q[mcl.df.q$FunGroup != "Native grass",], aes(x = FunGroup, y = c, col = FunGroup, shape = FunGroup)) +
  geom_boxplot(outliers = F) +
  geom_text(data = letters_df, aes(x = FunGroup, y = 0.26, label = label),
            vjust = 0, size = 6, col = "black") +
  geom_point(position = position_jitter(width = 0.35, height = 0), size = 1.8) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.position = "none",
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_blank(),
    
  ) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_y_continuous(limits = c(0.02, 0.28), breaks = c(0.1, 0.2)) +
  labs(x = "",
       y = "Probability of \n spatial dispersal")

fig1.c <- ggplot(mcl.df.q[mcl.df.q$FunGroup != "Native grass",], aes(x = FunGroup, y = s, col = FunGroup, shape = FunGroup)) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35, height = 0), size = 1.8) +
  geom_text(data = letters_df, aes(x = FunGroup, y = 0.85, label = label),
            vjust = 0, size = 6, col = "black") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.position = "none",
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 15),
  ) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_x_discrete(labels = c("Non-native\ngrass", "Non-native\nforb", "Native\nforb")) +
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8), limits = c(0.1,0.9)) +
  labs(x = "",
       y = "Probability of \n temporal dispersal")

fig1.bc <- ggarrange(fig1.b, fig1.c, ncol = 1, nrow = 2, labels = c("(b)", "(c)"), heights = c(1, 1.3), font.label = list(size = 17))
          
fig1 <- ggarrange(fig1.a, fig1.bc, ncol = 2, nrow = 1, labels = c("(a)", "", ""), widths = c(1.7, 1.3), font.label = list(size = 17)) + 
  bgcolor("white") +
  border(color = "white")

ggsave("Manuscript/Figures/Fig-1.pdf", fig1, height = 6, width = 10, units = "in", dpi = 600)

##### Models ####
par(mfrow = c(2,2))

cor.test(mcl.df.q$s, mcl.df.q$c, method = "spearman")
test <- lm(c ~ s, data = mcl.df.q)
plot(test)

Q1s.m <- lm(s ~ FunGroup, data = mcl.df.q[mcl.df.q$FunGroup != "Native grass",])
plot(Q1s.m)
pairs(emmeans(Q1s.m, ~ FunGroup), adjust = "BH")

Q1c.m <- lm(c ~ FunGroup, data = mcl.df.q[mcl.df.q$FunGroup != "Native grass",])
plot(Q1c.m)
emmeans(Q1c.m, ~ FunGroup)
pairs(emmeans(Q1c.m, ~ FunGroup), adjust = "BH")

par(mfrow = c(1,1))


#### Supplemental Table S1 ####

# dplyr::count(mcl.df.q, FunGroup)
# dplyr::count(mcl.df.q.trait, FunGroup)
# 
# acc <- read.csv("/Users/marina.laforgia/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/20230530_Seeds_All-Accessions.csv")
# ID <- read.csv("/Users/marina.laforgia/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/20230801_Seed-Traits_clean_ID.csv")
# 
# ID <- merge(ID[,c(2:5)], mcl.df.q.trait, by.x = c("Species", "site"), by.y = c("Species_Name", "site"), all.x = F, all.y = T)
# ID <- merge(acc[,c(1,14)], ID, by = "ID", all.x = F, all.y = T)
# ID$year <- as.numeric(ID$year)

#### Figure 2: Cor Plots ####

##### Spatial ####
colonization <- c("c", "new.mass",  "shape", "set.time.mpsec", "height.cm", "size.mm", "ldd.natural", "wing.loading")

trait.c <- data.frame(trait = character(), p.value = numeric(), cor = numeric(), FunGroup = character())

for(i in colonization[-1]){
  for(j in unique(mcl.df.q.trait[mcl.df.q.trait$FunGroup != "Native grass",]$FunGroup)) {
    
    tmp <- cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$c, mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i])
    
    trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = j, ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
    
  }
    tmp <- cor.test(mcl.df.q.trait$c, mcl.df.q.trait[,i])
    
  trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = "All", ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
}  

trait.c$trait <- factor(trait.c$trait, levels = trait.c[trait.c$FunGroup == "All",]$trait[order(abs(trait.c[trait.c$FunGroup == "All",]$cor))])

trait.c$FunGroup <- factor(trait.c$FunGroup, levels = c("Non-native grass", "Non-native forb", "Native forb", "All"))

c <- ggplot(trait.c, aes(x = cor, y = trait, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7), aes(fill = FunGroup)) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper), width = 0.1, position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(linewidth = 1, fill = NA)
  ) +
  scale_x_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  scale_shape_manual(values = c(15, 17, 19, 23)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  scale_y_discrete(labels = c("Mass", "Height", "Settling\n speed",  "Size", "Wing\n loading", "Shape",  "Dispersal\n potential")) + 
  labs(x = "Correlation", title = "Spatial dispersal")  


##### Temporal ####
survival <- c("s", "new.mass",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width")

trait.s <- data.frame(trait = character(), p.value = numeric(), cor = numeric(), FunGroup = character(), ci.lower = numeric(), ci.upper = numeric())

for(i in survival[-1]){
  for(j in unique(mcl.df.q.trait[mcl.df.q.trait$FunGroup != "Native grass",]$FunGroup)) {

    tmp <- cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$s, mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i])

    trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = j, ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))

  }
  tmp <- cor.test(mcl.df.q.trait$s, mcl.df.q.trait[,i])

  trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = "All", ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
}

trait.s$FunGroup <- factor(trait.s$FunGroup, levels = c("Non-native grass", "Non-native forb", "Native forb", "All"))

trait.s$trait <- factor(trait.s$trait, levels = trait.s[trait.s$FunGroup == "All",]$trait[order(abs(trait.s[trait.s$FunGroup == "All",]$cor))])
 
s <- ggplot(trait.s, aes(x = cor, y = trait, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7), aes(fill = FunGroup)) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper), width = 0.1, position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(linewidth = 1, fill = NA)
  ) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_shape_manual(values = c(15, 17, 19, 23)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  scale_y_discrete(labels = c("%C", "Mass", "Coat\n thickness", "Coat\n permeability", "Size", "%N", "Shape")) +
  labs(x = "Correlation", title = "Temporal dispersal")  

cor.plot <- ggarrange(s,c, widths = c(0.67,1), labels = c("(a)", "(b)"))
  
ggsave("Manuscript/Figures/Fig-2.pdf", cor.plot, height = 6, width = 9, units = "in", dpi = 600)


#### Figure 3: HMM & traits  ####
#inspect trait.c and trait.c for which traits

##### Shape ####
a <- ggplot(mcl.df.q.trait, aes(x = shape, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",], method = "lm", aes(col = FunGroup), se = F, linewidth = 1.3) +
  geom_smooth(inherit.aes = F, aes(x = shape, y = c), method = "lm", col = "black", se = F, linewidth = 1.3) +  
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.038,0.28)) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(y = "Spatial dispersal", x = "Shape")

f <- ggplot(mcl.df.q.trait, aes(x = shape, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",], method = "lm", aes(col = FunGroup), se = F, linewidth = 1.3) + 
  geom_smooth(inherit.aes = F, aes(x = shape, y = s), method = "lm", col = "black", se = F, size = 1.3) +  
  geom_point(size = 2, aes(fill = FunGroup)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(y = "Temporal dispersal", x = "Shape")

##### Size ####

b <- ggplot(mcl.df.q.trait, aes(x = size.mm, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",], method = "lm", aes(col = FunGroup), linetype = 2, se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = size.mm, y = s), method = "lm", col = "black", se = F, linetype = 2, size = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "Log size (mm)")

g <- ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "Non-native grass",], method = "lm", aes(col = FunGroup), se = F, linetype = 2, linewidth = 1.3) + 
  geom_smooth(inherit.aes = F, aes(x = size.mm, y = c), method = "lm", col = "black", se = F, size = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),    
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.038,0.28)) +
  scale_linetype_manual(values = c(1,2)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "Log size (mm)")

##### Coat thick ####

c <- ggplot(mcl.df.q.trait, aes(x = coat.thick.per.width, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "Native forb",], aes(col = FunGroup, linetype = FunGroup), method = "lm", se = F, linewidth = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_linetype_manual(values = c(1,2)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "Seed coat ratio")

  
##### N ####

d <- ggplot(mcl.df.q.trait, aes(x = prop.N, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "Native forb",], method = "lm", aes(fill = FunGroup), se = F, linewidth = 1.3) + 
  geom_smooth(inherit.aes = F, aes(x = prop.N, y = s), method = "lm", col = "black", se = F, size = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "%N")

##### Mass ####
i <- ggplot(mcl.df.q.trait, aes(x = new.mass, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native grass",], method = "lm", col = "#1b9e77", se = F, linewidth = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(y = "Temporal dispersal", x = "Log mass (mg)")

##### Disp ####

e <- ggplot(mcl.df.q.trait, aes(x = ldd.natural, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",], method = "lm", col = "#7570b3", se = F, linewidth = 1.3) +
  geom_smooth(inherit.aes = F, aes(x = ldd.natural, y = c), method = "lm", col = "black", se = F, linewidth = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(color = "white")
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.038,0.28)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "Dispersal potential", y = "Spatial dispersal")


##### Wing loading ####

o <- ggplot(mcl.df.q.trait, aes(x = wing.loading, y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",], method = "lm", fill = "#7570b3", se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = wing.loading, y = c), method = "lm", col = "black", se = F, linewidth = 1.3) +
  geom_point(size = 2, aes(fill = FunGroup)) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "Log wing loading")


##### Panel Plot ####
leg <- get_legend(
  ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
    theme_classic() +
    geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2, se = F) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    ) +
    geom_point(size = 3, aes(fill = FunGroup)) +
    scale_color_manual(values = fungroup_cols) +
    scale_shape_manual(values = fungroup_shapes) +
    scale_fill_manual(values = fungroup_cols)
  )

trait.panel.s <- ggarrange(e, a, o, g, leg, ncol = 5, nrow = 1, 
                           labels = c("(a)", "(b)", "(c)" , "(d)"), 
                           widths = c(1.07, 0.85, 0.85, 0.85, 0.85), 
                           heights = 1, 
                           label.x = c(0.22,0.04,0.04,0.04), 
                           label.y = 0.96) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.s <- annotate_figure(trait.panel.s, top = text_grob("Spatial dispersal", color = "black", face = "bold", size = 18)) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.t <- ggarrange(f,d,b,c, i, ncol = 5, nrow = 1, 
                           labels = c("(e)", "(f)", "(g)", "(h)", "(i)"), 
                           widths = c(1.07, 0.85, 0.85, 0.85, 0.85), 
                           heights = 1, 
                           label.x = c(0.22,0.04,0.03,0.03,0.04), 
                           label.y = 0.96) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.t <- annotate_figure(trait.panel.t, top = text_grob("Temporal dispersal", color = "black", face = "bold", size = 18)) + 
  bgcolor("white") +
  border(color = "white")

trait.panel <- ggarrange(trait.panel.s, trait.panel.t, ncol = 1, nrow = 2, align = "hv") #, heights = c(0.9, 1))
                           
ggsave("Manuscript/Funct-Ecol/Revision/Final-Figs-PDF-600-DPI/Fig-3.pdf", trait.panel, height = 6, width = 12, units = "in", dpi = 600)

#### Linear models ####
par(mfrow = c(2,2))

##### Shape ####
m.shape.s <- lm(s ~ shape, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",])
plot(m.shape.s)
summary(m.shape.s)
  
m.shape.c <- lm(c ~ shape, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",])
plot(m.shape.c)
summary(m.shape.c)

m.shape.s <- lm(s ~ shape, data = mcl.df.q.trait)
plot(m.shape.s)
summary(m.shape.s)

m.shape.c <- lm(c ~ shape, data = mcl.df.q.trait)
plot(m.shape.c)
summary(m.shape.c)

##### Size ####
m.size.s <- lm(s ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb", ])
plot(m.size.s)
summary(m.size.s)
  
m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",])
plot(m.size.c)
summary(m.size.c)

m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native forb",])
plot(m.size.c)
summary(m.size.c)

m.size.s <- lm(s ~ size.mm, data = mcl.df.q.trait)
plot(m.size.s)
summary(m.size.s)

m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait)
plot(m.size.c)
summary(m.size.c)


##### Wing loading ####
m.wing.c <- lm(c ~ wing.loading, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",])
plot(m.wing.c)
summary(m.wing.c)

m.wing.c <- lm(c ~ wing.loading, data = mcl.df.q.trait)
plot(m.wing.c)
summary(m.wing.c)

##### Coat ####
m.coat.s <- lm(s ~ coat.thick.per.width, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native forb",])
plot(m.coat.s)
summary(m.coat.s)

m.coat.s.EG <- lm(s ~ coat.thick.per.width, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native grass",])
plot(m.coat.s.EG)
summary(m.coat.s.EG) # marginal relationship

##### N ####
m.N.c.EF <- lm(s ~ prop.N, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native forb",])
plot(m.N.c.EF)
summary(m.N.c.EF)

m.N.c.EG <- lm(s ~ prop.N, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native grass",])
plot(m.N.c.EG)
summary(m.N.c.EG) # GASVEN and BROMAD have high leverage

m.N.c.EF <- lm(s ~ prop.N, data = mcl.df.q.trait)
plot(m.N.c.EF)
summary(m.N.c.EF)

##### Disp ####
m.disp.c <- lm(c ~ ldd.natural, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Native forb",])
plot(m.disp.c)
summary(m.disp.c)

m.disp.c <- lm(c ~ ldd.natural, data = mcl.df.q.trait)
plot(m.disp.c)
summary(m.disp.c)

##### Mass ####
m.mass.s <- lm(s ~ new.mass, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "Non-native grass",])

plot(m.mass.s)
summary(m.mass.s)

par(mfrow = c(1,1))

#### Phylogenetic analyses ####

##### Prep data ####

# lots of reconciling synonyms here

mcl.df.q.all$Species_Name_syn <- mcl.df.q.all$Species_Name
levels(mcl.df.q.all$Species_Name_syn) <- c(levels(mcl.df.q.all$Species_Name_syn), 
                                           "Phlox gracilis",
                                           "Lysimachia arvensis",
                                           "Stylocline filaginea",
                                           "Sida diploscypha",
                                           "Lepidostephanus madioides",
                                           "Centaurium trichanthum",
                                           "Microseris lindleyi")

# synonyms
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Microsteris gracilis",]$Species_Name_syn <- "Phlox gracilis"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Anagallis arvensis",]$Species_Name_syn <- "Lysimachia arvensis"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Ancistrocarphus filagineus",]$Species_Name_syn <- "Stylocline filaginea"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Sidalcea diploscypha",]$Species_Name_syn <- "Sida diploscypha"
#mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Zeltnera trichantha",]$Species_Name_syn <- "Centaurium trichanthum"
#mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Uropappus lindleyi",]$Species_Name_syn <- "Microseris lindleyi"

# missing families
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Acmispon brachycarpus",]$family <- "Fabaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Cerastium glomeratum",]$family <- "Caryophyllaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Castilleja rubicundula",]$family <- "Orobanchaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Navarretia jepsonii",]$family <- "Polemoniaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Leptosiphon bicolor",]$family <- "Polemoniaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Lactuca serriola",]$family <- "Asteraceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Croton setiger",]$family <- "Euphorbiaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Triphysaria eriantha",]$family <- "Orobanchaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Eriogonum vimineum",]$family <- "Polygonaceae"
mcl.df.q.all[mcl.df.q.all$Species_Name_syn == "Galium aparine",]$family <- "Rubiaceae"

mcl.df.q.all$species <- gsub(" ", "_", mcl.df.q.all$Species_Name_syn)

##### Build tree ####
gen.list <- read.csv("Data/Tree/plant_genus_list.csv")
tree <- read.tree("Data/Tree/plant_megatree.tre")


result <- phylo.maker(mcl.df.q.all$species, tree = tree, gen.list = gen.list)

tree2 <- result$phylo

is.ultrametric(tree2)
max(branching.times(tree2)) 

##### Fig S1: Tree ####

# Align species info to tree
df <- mcl.df.q.all[match(tree2$tip.label, mcl.df.q.all$species), ]

# Define functional group colors
fg_cols <- c(
  "non-native forb"  = "#d95f02",
  "non-native grass" = "#1b9e77",
  "native forb"      = "#7570b3",
  "native grass"     = "darkgreen"
)

# Assign a color to each tip based on its FunGroup
tip_cols <- fg_cols[df$FunGroup]

# Make sure data frame has a "family" column aligned with the tree
fam <- df$family
fam_levels <- unique(fam)

tree2$tip.label[tree2$tip.label == "Sida_diploscypha"] <- "Sidalcea_diploscypha"

pdf("Manuscript/Figures/Fig-S1.pdf", width = 8, height = 10)

plot(
  tree2,
  cex = 0.7,
  tip.color = tip_cols,
  edge.width = 1,
  label.offset = 0.3
)
axisPhylo()
mtext("Time (millions of years ago)", side = 1, line = 2.5, cex = 0.9)

# Add legend
legend(
  "topleft",
  legend = names(fg_cols),
  col = fg_cols,
  pch = 19,
  cex = 0.8,
  bty = "n",
  y.intersp = 0.9 
)

for (f in fam_levels) {
  tips <- tree2$tip.label[fam == f]
  if (length(tips) > 1) {
    node <- getMRCA(tree2, tips)
    if (!is.null(node) && !is.na(node)) {
      nodelabels(
        text = f,
        node = node,
        frame = "none",
        cex = 0.7,
        col = "gray30",
        adj = c(1.2, -0.3)
      )
    }
  }
}

singletons <- fam_levels[sapply(fam_levels, function(f)
  sum(fam == f) == 1)]

if (length(singletons) > 0) {
  for (f in singletons) {
    tip <- which(fam == f)
    tiplabels(
      text = f,
      tip = tip,
      frame = "none",
      cex = 0.7,
      col = "gray30",
      adj = c(2, -0.3)
    )
  }
}

dev.off()

##### Inspect Tree output ####
tree2$tip.label[tree2$tip.label == "Sidalcea_diploscypha"] <- "Sida_diploscypha"
tree.sp <- result$sp.list # most at species, a lot at genus, one at family

# make sure names match
length(tree2$tip.label)                    # total in tree: 60
nrow(mcl.df.q.all)                        # total in data: 60
length(intersect(tree2$tip.label, mcl.df.q.all$species))  #: 60


##### PGLS c-s ####
# now lets look a phylogenetic signal
# colnames(mcl.df.q.all)[55] <- "species"
# mcl.df.q.all$species <- gsub(" ", "_", mcl.df.q.all$species)
tree2$node.label <- NULL

## no phylogenetic signal in C or S
phylosig(tree2, mcl.df.q.all$c, method = "lambda", test = TRUE)
phylosig(tree2, mcl.df.q.all$s, method = "lambda", test = TRUE)

sum(mcl.df.q.all$species %in% tree2$tip.label)           # should be 60
identical(sort(mcl.df.q.all$species), sort(tree2$tip.label))

mcl.df.q.all <- subset(mcl.df.q.all, !is.na(c) & !is.na(s))
common <- intersect(tree2$tip.label, mcl.df.q.all$species)
tree2 <- drop.tip(tree2, setdiff(tree2$tip.label, common))

comp <- comparative.data(tree2, mcl.df.q.all, names.col = "species", na.omit = F) 
pgls_fit <- pgls(c ~ s, data = comp, lambda = "ML") # lambda based on maximum likelihood
summary(pgls_fit)


##### Phylosigs for traits ####
phylosig(tree2, mcl.df.q.all$shape, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$new.mass, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$coat.thick.per.width, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$size.mm, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$prop.N, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$ldd.natural, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$ldd.natural, method = "lambda", test = TRUE) # yes
phylosig(tree2, mcl.df.q.all$wing.loading, method = "lambda", test = TRUE) # yes

traits <- c("shape", "new.mass", "coat.thick.per.width", "size.mm", "prop.N", "ldd.natural", "wing.loading", "prop.C", "coat.perm.perc", "set.time.mpsec", "height.cm")

# initialize results dataframe
results <- data.frame(
  trait = character(),
  lambda = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (tr in traits) {
  
  # extract and name trait vector
  x <- mcl.df.q.all[[tr]]
  names(x) <- mcl.df.q.all$species
  
  # remove NAs
  x <- x[!is.na(x)]
  
  # match species with tree
  common <- intersect(tree2$tip.label, names(x))
  x <- x[common]
  tree.pruned <- drop.tip(tree2, setdiff(tree2$tip.label, common))
  
  # run phylosig
  ps <- phylosig(
    tree.pruned,
    x,
    method = "lambda",
    test = TRUE
  )
  
  # store results
  results <- rbind(
    results,
    data.frame(
      trait = tr,
      lambda = round(ps$lambda, 3),
      p_value = round(ps$P, 3),
      stringsAsFactors = FALSE
    )
  )
}

results <- results[order(results$lambda), ]

#### Fig S2: Correlation ####
mcl.sim <- readRDS("Intermediates/McL-sim-154-correlation.RDS")

mcl.sim <- as.data.frame(mcl.sim)

mcl.sim <- filter(mcl.sim, iter < 150)

sim.df <- data.frame(cor = NA, p.value = NA)

for(i in 1:50000){
  tmp <- slice_sample(mcl.sim, n = 60, replace = F)
  tmp.test <- cor.test(tmp$s, tmp$c)
  sim.df[i,]$cor <- tmp.test$estimate
  sim.df[i,]$p.value <- tmp.test$p.value
}

figs1 <- ggplot(sim.df, aes(x = cor)) +
  geom_histogram(bins = 20, col = "black", fill = "white") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  geom_vline(xintercept = -0.61, col = "red") +
  labs(x = "simulated correlation coefficient")

mean(sim.df$cor)
nrow(sim.df[sim.df$cor <= -0.61,])/nrow(sim.df) #p value < 0.0001

ggsave("Manuscript/Figures/Final-Figs-PDF-600-DPI/Fig-S2.pdf", figs1, height = 5, width = 7, units = "in", dpi = 600)

#### Fig S3: Species validation ####
mcl.sim.sp <- readRDS("/Users/marina.laforgia/Documents/USDA-PostDoc/Projects/Seed-Traits/Scripts/HMMs/McLaughlin/McL-sim-species-validation-50.RDS")

mcl.sim.sp <- as.data.frame(mcl.sim.sp)

mcl.sim.sp <- mcl.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(mcl.sim.sp)[4:9] <- paste(colnames(mcl.sim.sp)[4:9], "est", sep=".")

mcl.sim.sp <- merge(mcl.sim.sp, mcl.df.q, by = "Species_Name")

mcl.rmse <- mcl.sim.sp %>%
  #group_by(sim, patches) %>%
  summarize(s.bias = mean(abs(s - mean(s.est))),
            s.rmse = rmse(s, s.est),
            s.var = var(s, s.est),
            c.bias = mean(abs(c - mean(c.est))),
            c.rmse = rmse(c, c.est),
            c.var = var(c, c.est),
            g.bias = mean(abs(g - mean(g.est))),
            g.rmse = rmse(g, g.est),
            g.var = var(g, g.est)
            )

mcl.sim.sum <- mcl.sim.sp %>%
  dplyr::group_by(Species_Name) %>%
  dplyr::summarize(across(p0.est:s, ~mean(.x)))

a <- ggplot(mcl.sim.sum, aes(y = s.est, x = s)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point(size = 0.9) +
  theme_classic() +
    theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  scale_y_continuous(breaks = c(0.3, 0.6)) +
  scale_x_continuous(breaks = c(0.3, 0.6)) +
  labs(x = "real s", y = "mean simulated s")


b <- ggplot(mcl.sim.sum, aes(y = c.est, x = c)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point(size = 0.9) +
  theme_classic() +
    theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real c", y = "mean simulated c")

c <- ggplot(mcl.sim.sum, aes(y = g.est, x = g)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point(size = 0.9) +
  theme_classic() +
    theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real g", y = "mean simulated g")

fig.s3 <- ggarrange(a,b,c, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"))

ggsave("Manuscript/Figures/Final-Figs-PDF-600-DPI/Fig-S3.pdf", fig.s3, height = 2.5, width = 8, units = "in", dpi = 600)


#### Fig S4: Other traits ####

##### Mass ####

h <- ggplot(mcl.df.q.trait, aes(x = new.mass, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(y = "spatial dispersal", x = "log mass (mg)")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$new.mass)

##### Set speed ####
j <- ggplot(mcl.df.q.trait, aes(x = set.time.mpsec, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "settling speed (m/s)")

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$set.time.mpsec) #ns

##### C ####
k <- ggplot(mcl.df.q.trait, aes(x = prop.C, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "%C", y = "temporal dispersal")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$prop.C) #ns

##### Coat perm ####
l <- ggplot(mcl.df.q.trait, aes(x = coat.perm.perc, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "log coat permeability")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$coat.perm.perc) #ns

##### Height ####
m <- ggplot(mcl.df.q.trait, aes(x = height.cm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_point(size = 2, aes(fill = FunGroup)) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_color_manual(values = fungroup_cols) +
  scale_shape_manual(values = fungroup_shapes) +
  scale_fill_manual(values = fungroup_cols) +
  labs(x = "height (cm)")

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$height.cm) #ns

##### Panel Plot ####
leg <- get_legend(
  ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
    theme_classic() + 
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    geom_point(size = 2, aes(fill = FunGroup)) + 
    scale_color_manual(values = fungroup_cols) +
    scale_shape_manual(values = fungroup_shapes) +
    scale_fill_manual(values = fungroup_cols)
  )

s4.appen.panel <- ggarrange(h,j,m,k,l, leg, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)"), widths = c(1, 0.85, 0.85), heights = c(1,1,1), label.x = c(0.2,0.04,0.04,0.2,0.04,0.04), label.y = 0.96) +
  bgcolor("white") +
  border(color = "white")

ggsave("Manuscript/Figures/Final-Figs-PDF-600-DPI/Fig-S4.pdf", s4.appen.panel, height = 6, width = 8.5, units = "in", dpi = 600)


#### Fig S5: Trait correlations ####
traits <- c("shape", "size.mm", "wing.loading", "ldd.natural","new.mass", "prop.N", "coat.thick.per.width", "set.time.mpsec", "height.cm", "prop.C",  "coat.perm.perc")

      
corr <- round(cor(mcl.df.q.trait[, traits]), 1)
p.mat <- cor_pmat(mcl.df.q.trait[, traits])

fig.s5 <- ggcorrplot(corr, 
           p.mat = p.mat,
           #method = "circle", 
           hc.order = TRUE, 
           #outline.col = "white", 
           #type = "lower",
           lab = T, 
           lab_size = 3,
           insig = "blank",
           colors = c("red", "white", "blue"),
           show.diag = F) +
  scale_y_discrete(labels = c("coat perm", "mass", "wing load", "set speed",  "height", "disp potential", "shape", "size", "coat thick", "%N",  "%C")) +
  scale_x_discrete(labels = c("coat perm", "mass", "wing load", "set speed",  "height", "disp potential", "shape", "size", "coat thick", "%N",  "%C")) 

ggsave("Manuscript/Figures/Final-Figs-PDF-600-DPI/Fig-S5.pdf", fig.s5, dpi = 600, units = "in", height = 5, width = 5)
