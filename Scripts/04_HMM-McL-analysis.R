#### HMMs McLaughlin Manuscript ####

rm(list = ls())

#### Load libraries ####
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(Metrics)
library(ggcorrplot)

#### Read in data ####
mcl.df.q <- readRDS("Data/McL-species-boot.RDS") # quadrat level HMM
metadata <- read.csv("Data/HMM-meta-mcl.csv")
seed.traits <- read.csv("Data/McL_seed-traits.csv")

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

mcl.df.q <- filter(mcl.df.q, FunGroup != "Native Grass") # one species lost

mcl.df.q$FunGroup <- recode_factor(mcl.df.q$FunGroup, 'Exotic Forb' = "non-native forb", 'Exotic Grass' = "non-native grass", 'Native Forb' = "native forb")

colnames(mcl.df.q)[2:6] <- c("p0", "g", "c", "s", "r")

#### Seed trait prep ####

# Normalize seed.traits
#hist(log(seed.traits$wing.loading))
seed.traits$wing.loading <- log(seed.traits$wing.loading)

#hist(log(seed.traits$coat.perm.perc))
seed.traits$coat.perm.perc <- log(seed.traits$coat.perm.perc)

#hist(log(seed.traits$morph.mass.mg))
seed.traits$morph.mass.mg <- log(seed.traits$morph.mass.mg)

#hist(log(seed.traits$chem.mass.mg))
seed.traits$chem.mass.mg <- log(seed.traits$chem.mass.mg)

# hist(log(seed.traits$new.mass))
seed.traits$new.mass <- log(seed.traits$new.mass)

#hist(log(seed.traits$size.mm))
seed.traits$size.mm <- log(seed.traits$size.mm)

mcl.df.q.trait <- merge(mcl.df.q, seed.traits, by.x = c("Species_Name", "FunGroup", "new.code"), by.y = c("Species", "FunGroup", "new.code"), all = F)

#### Figure 1: Trade-off ####
##### Figure ####

# prep figure to graph CESO datapoint on top (makes it easier to see with the seed pic)
fig1.df <- mcl.df.q
fig1.df$s1 <- NA
fig1.df$c1 <- NA
fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$s1 <- fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$s
fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$c1 <- fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$c

fig1 <- ggplot(fig1.df, aes(x = s, y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(inherit.aes = F, aes(x = s, y = c), method = "lm", col = "black", se = F) + 
  geom_point(size = 3) +
  #geom_text(aes(label = new.code)) +
  geom_point(aes(x = s1, y = c1, col = FunGroup, shape = FunGroup), size = 3) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.position = "right",
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
  ) +
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "probability of temporal dispersal",
       y = "probability of spatial dispersal")

ggsave("Manuscript/Figures/C-S.tiff", fig1, height = 5, width = 7, units = "in", dpi = 600)

##### Models ####
par(mfrow = c(2,2))

cor.test(mcl.df.q$s, mcl.df.q$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "non-native forb",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "non-native forb",]$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "native forb",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "native forb",]$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "non-native grass",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "non-native grass",]$c)

Q1s.m <- lm(s ~ FunGroup, data = mcl.df.q)
plot(Q1s.m)
emmeans(Q1s.m, ~ FunGroup)
pairs(emmeans(Q1s.m, ~ FunGroup), adjust = "BH")

Q1c.m <- lm(c ~ FunGroup, data = mcl.df.q)
plot(Q1c.m)
emmeans(Q1c.m, ~ FunGroup)
pairs(emmeans(Q1c.m, ~ FunGroup), adjust = "BH")

par(mfrow = c(1,1))

#### Figure 2: Cor Plots ####

##### Spatial ####
colonization <- c("c", "new.mass",  "shape", "set.time.mpsec", "height.cm", "size.mm", "ldd.natural", "wing.loading")

trait.c <- data.frame(trait = character(), p.value = numeric(), cor = numeric(), FunGroup = character())

for(i in colonization[-1]){
  for(j in unique(mcl.df.q.trait$FunGroup)) {
    
    tmp <- cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$c, mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i])
    
    trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = j, ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
    
  }
    tmp <- cor.test(mcl.df.q.trait$c, mcl.df.q.trait[,i])
    
  trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = "all", ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
}  

trait.c$trait <- factor(trait.c$trait, levels = trait.c[trait.c$FunGroup == "all",]$trait[order(abs(trait.c[trait.c$FunGroup == "all",]$cor))])

trait.c$FunGroup <- factor(trait.c$FunGroup, levels = c("non-native grass", "non-native forb", "native forb", "all"))

c <- ggplot(trait.c, aes(x = cor, y = trait, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_rect(mapping = aes(xmin = -1, xmax = 1, 
  #                         ymin = 2.57, ymax = 3.45), col = "grey86", fill =    #"grey86") +
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
  scale_y_discrete(labels = c("mass", "height", "settling\n speed",  "size", "wing\n loading", "shape",  "dispersal\n potential")) + 
  labs(x = "correlation", title = "Spatial Dispersal")  


##### Temporal ####
survival <- c("s", "new.mass",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width")

trait.s <- data.frame(trait = character(), p.value = numeric(), cor = numeric(), FunGroup = character(), ci.lower = numeric(), ci.upper = numeric())

for(i in survival[-1]){
  for(j in unique(mcl.df.q.trait$FunGroup)) {
      
    tmp <- cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$s, mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i])
    
    trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = j, ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
    
  }
  tmp <- cor.test(mcl.df.q.trait$s, mcl.df.q.trait[,i])
    
  trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = "all", ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
}  

trait.s$FunGroup <- factor(trait.s$FunGroup, levels = c("non-native grass", "non-native forb", "native forb", "all"))

trait.s$trait <- factor(trait.s$trait, levels = trait.s[trait.s$FunGroup == "all",]$trait[order(abs(trait.s[trait.s$FunGroup == "all",]$cor))])
 
s <- ggplot(trait.s, aes(x = cor, y = trait, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  # geom_rect(mapping = aes(xmin = -1, xmax = 1, 
  #                         ymin = 3.55, ymax = 4.45), col = "grey86", fill = "grey86") +
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
  scale_y_discrete(labels = c("%C", "mass", "coat\n thickness", "coat\n permeability", "size", "%N", "shape")) +
  labs(x = "correlation", title = "Temporal Dispersal")  

cor.plot <- ggarrange(s,c, widths = c(0.67,1), labels = c("(a)", "(b)"))
  
ggsave("cor-plot.png", cor.plot, height = 6, width = 9, units = "in", dpi = 600)


#### Figure 3: HMM & traits  ####

##### Shape ####
 
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$shape)  #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$shape) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$shape) #ns

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$shape) #sig

a <- ggplot(mcl.df.q.trait, aes(x = shape, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", aes(col = FunGroup), se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = shape, y = c), method = "lm", col = "black", se = F, linewidth = 1.3) +  
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$shape) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$shape) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$shape) #ns

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$shape) # sig

f <- ggplot(mcl.df.q.trait, aes(x = shape, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", aes(col = FunGroup), se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = shape, y = s), method = "lm", col = "black", se = F, size = 1.3) +  
  #geom_text(aes(label = new.code.x)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "temporal dispersal")

##### Size ####
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$size.mm) #marg

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$size.mm) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$size.mm) #ns

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$size.mm) # marg

b <- ggplot(mcl.df.q.trait, aes(x = size.mm, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", aes(col = FunGroup), linetype = 2, se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = size.mm, y = s), method = "lm", col = "black", se = F, linetype = 2, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  #scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "log size (mm)")
 
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$size.mm) #marg

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$size.mm) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$size.mm) #marg

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$size.mm) # sig

g <- ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "non-native grass",], method = "lm", aes(col = FunGroup), se = F, linetype = 2, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = size.mm, y = c), method = "lm", col = "black", se = F, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#d95f02", "#7570b3","#1b9e77")) +
  scale_shape_manual(values = c( 17, 19, 15)) +
  labs(x = "log size (mm)")

##### Coat thick ####
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$coat.thick.per.width) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$coat.thick.per.width) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$coat.thick.per.width) # marg

# cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$prop.N, 
#          mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$coat.thick.per.width)

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$coat.thick.per.width) #ns

c <- ggplot(mcl.df.q.trait, aes(x = coat.thick.per.width, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "native forb",], aes(col = FunGroup, linetype = FunGroup), method = "lm", se = F, linewidth = 1.3) +
  #geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",], method = "lm", col = "#1b9e77", se = F,  linetype = 2) + 
  geom_point(size = 2) +
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
  #scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_linetype_manual(values = c(1,2)) +
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "seed coat ratio")

  
##### N ####

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$prop.N) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$prop.N) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$prop.N) #sig

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$prop.N) #sig

d <- ggplot(mcl.df.q.trait, aes(x = prop.N, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "native forb",], method = "lm", aes(fill = FunGroup), se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = prop.N, y = s), method = "lm", col = "black", se = F, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  #scale_fill_manual(values = c("#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "%N")



##### Mass ####
i <- ggplot(mcl.df.q.trait, aes(x = new.mass, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",], method = "lm", col = "#1b9e77", se = F, linewidth = 1.3) +
  geom_point(size = 2) +
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
  #scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.16,0.85)) +
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  scale_shape_manual(values = c(15, 17, 19)) +
  labs(y = "temporal dispersal", x = "log mass (mg)")

##### Disp ####
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$ldd.natural) #ns


cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$ldd.natural) #sig

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$ldd.natural) #sig

e <- ggplot(mcl.df.q.trait, aes(x = ldd.natural, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", col = "#7570b3", se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = ldd.natural, y = c), method = "lm", col = "black", se = F, linewidth = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "dispersal potential", y = "spatial dispersal")


##### Wing loading ####
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$wing.loading) #ns


cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$wing.loading) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$wing.loading) #ns

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$wing.loading) #sig

o <- ggplot(mcl.df.q.trait, aes(x = wing.loading, y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", fill = "#7570b3", se = F, linewidth = 1.3) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = wing.loading, y = c), method = "lm", col = "black", se = F, linewidth = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "log wing loading")


##### Panel Plot ####
leg <- get_legend(
  ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
    theme_classic() +
    geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2, se = F) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    ) +
    geom_point(size = 3) +
    scale_color_manual(values = c("#d95f02", "#1b9e77","#7570b3")) +
    scale_shape_manual(values = c(17, 15, 19))
  )
# 
# trait.panel <- ggarrange(a,g,e,o,f,b,c,d, ncol = 4, nrow = 2, labels = c("(a)", "(b)", "(c)" , "(d)", "(e)", "(f)", "(g)", "(h)"), widths = c(1, 0.85, 0.85, 0.85), heights = c(1,1), label.x = c(0.2,0.04,0.04,0.04,0.2,0.04,0.04,0.04), label.y = 0.96, common.legend = TRUE,legend="bottom")

trait.panel.s <- ggarrange(e, a, o, g, leg, ncol = 5, nrow = 1, labels = c("(a)", "(b)", "(c)" , "(d)"), widths = c(1.07, 0.85, 0.85, 0.85, 0.85), heights = 1, label.x = c(0.2,0.04,0.04,0.04), label.y = 0.96) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.s <- annotate_figure(trait.panel.s, top = text_grob("Spatial Dispersal", 
               color = "black", face = "bold", size = 18)) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.t <- ggarrange(f,d,b,c, i, ncol = 5, nrow = 1, labels = c("(e)", "(f)", "(g)", "(h)", "(i)"), widths = c(1.07, 0.85, 0.85, 0.85, 0.85), heights = 1, label.x = c(0.2,0.04,0.04,0.04,0.04), label.y = 0.96) + #, common.legend = TRUE, legend = "bottom"
  bgcolor("white") +
  border(color = "white")

trait.panel.t <- annotate_figure(trait.panel.t, top = text_grob("Temporal Dispersal", color = "black", face = "bold", size = 18)) + 
  bgcolor("white") +
  border(color = "white")

trait.panel <- ggarrange(trait.panel.s, trait.panel.t, ncol = 1, nrow = 2, align = "hv") #, heights = c(0.9, 1))
                           
ggsave("Manuscript/Figures/Fig-3.png", trait.panel, height = 6, width = 12, units = "in", dpi = 600)

#### Linear models ####
par(mfrow = c(2,2))

##### Shape ####
m.shape.s <- lm(s ~ shape, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.shape.s)
summary(m.shape.s)
  
m.shape.c <- lm(c ~ shape, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.shape.c)
summary(m.shape.c)

m.shape.s <- lm(s ~ shape, data = mcl.df.q.trait)
plot(m.shape.s)
summary(m.shape.s)

m.shape.c <- lm(c ~ shape, data = mcl.df.q.trait)
plot(m.shape.c)
summary(m.shape.c)

##### Size ####
m.size.s <- lm(s ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.size.s)
summary(m.size.s)
  
m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.size.c)
summary(m.size.c)

m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",])
plot(m.size.c)
summary(m.size.c)

m.size.s <- lm(s ~ size.mm, data = mcl.df.q.trait)
plot(m.size.s)
summary(m.size.s)

m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait)
plot(m.size.c)
summary(m.size.c)


##### Wing loading ####
m.wing.c <- lm(c ~ wing.loading, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.wing.c)
summary(m.wing.c)

m.wing.c <- lm(c ~ wing.loading, data = mcl.df.q.trait)
plot(m.wing.c)
summary(m.wing.c)

##### Coat ####
m.coat.s <- lm(s ~ coat.thick.per.width, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",])
plot(m.coat.s)
summary(m.coat.s)

m.coat.s.EG <- lm(s ~ coat.thick.per.width, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",])
plot(m.coat.s.EG)
summary(m.coat.s.EG)

##### N ####
m.N.c.EF <- lm(s ~ prop.N, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",])
plot(m.N.c.EF)
summary(m.N.c.EF)

m.N.c.EG <- lm(s ~ prop.N, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",])
plot(m.N.c.EG)
summary(m.N.c.EG)

m.N.c.EF <- lm(s ~ prop.N, data = mcl.df.q.trait)
plot(m.N.c.EF)
summary(m.N.c.EF)

##### Disp ####
m.disp.c <- lm(c ~ ldd.natural, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.disp.c)
summary(m.disp.c)

m.disp.c <- lm(c ~ ldd.natural, data = mcl.df.q.trait)
plot(m.disp.c)
summary(m.disp.c)

##### Mass ####
m.mass.s <- lm(s ~ new.mass, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",])
plot(m.mass.s)
summary(m.mass.s)

par(mfrow = c(1,1))

#### Fig S1: Correlation ####
mcl.sim <- readRDS("Data/McL-sim-154-correlation.RDS")

mcl.sim <- as.data.frame(mcl.sim)

mcl.sim <- filter(mcl.sim, iter < 150)

sim.df <- data.frame(cor = NA, p.value = NA)

for(i in 1:50000){
  tmp <- slice_sample(mcl.sim, n = 59, replace = F)
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
  geom_vline(xintercept = -0.536, col = "red") +
  labs(x = "simulated correlation coefficient")

mean(sim.df$cor)
nrow(sim.df[sim.df$cor <= -0.536,])/nrow(sim.df) #p value < 0.0001

ggsave("Manuscript/Figures/figs1.png", figs1, height = 5, width = 7, units = "in", dpi = 600)

#### Fig S2: Species validation ####
mcl.sim.sp <- readRDS("Data/McL-sim-species-validation-50.RDS")

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

fig.s2 <- ggarrange(a,b,c, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"))

ggsave("Manuscript/Figures/FigS2-species-validation.png", fig.s2, height = 2.5, width = 8, units = "in", dpi = 600)


#### Fig S3: Other traits ####

##### Mass ####

h <- ggplot(mcl.df.q.trait, aes(x = new.mass, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(y = "spatial dispersal", x = "log mass (mg)")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$new.mass)

##### Set speed ####
j <- ggplot(mcl.df.q.trait, aes(x = set.time.mpsec, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "settling speed (m/s)")

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$set.time.mpsec) #ns

##### C ####
k <- ggplot(mcl.df.q.trait, aes(x = prop.C, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "%C", y = "temporal dispersal")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$prop.C) #ns

##### Coat perm ####
l <- ggplot(mcl.df.q.trait, aes(x = coat.perm.perc, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "log coat permeability")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$coat.perm.perc) #ns

##### Height ####
m <- ggplot(mcl.df.q.trait, aes(x = height.cm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
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
   scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
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
    geom_point(size = 2) + 
   scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19))
  )

s3.appen.panel <- ggarrange(h,j,m,k,l, leg, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)"), widths = c(1, 0.85, 0.85), heights = c(1,1,1), label.x = c(0.2,0.04,0.04,0.2,0.04,0.04), label.y = 0.96) +
  bgcolor("white") +
  border(color = "white")

ggsave("Manuscript/Figures/Fig-s3.png", s3.appen.panel, height = 6, width = 8.5, units = "in", dpi = 600)


#### Fig S4: Trait correlations ####
traits <- c("shape", "size.mm", "wing.loading", "ldd.natural","new.mass", "prop.N", "coat.thick.per.width", "set.time.mpsec", "height.cm", "prop.C",  "coat.perm.perc")

      
corr <- round(cor(mcl.df.q.trait[, traits]), 1)
p.mat <- cor_pmat(mcl.df.q.trait[, traits])

fig.s4 <- ggcorrplot(corr, 
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

ggsave("Manuscript/Figures/Fig-s4.png", fig.s4, dpi = 600, units = "in", height = 5, width = 5)
