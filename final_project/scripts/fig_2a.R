# fig_2a.R
#- Fig2A: Superclusters enrichment significance of traits' SNP-heritability

library(tidyverse)
library(ggplot2)
library(here)
library(readxl)
library(scico)
source("../scripts/order_data.R"))

dat1 <- read_xlsx("../data/Supplementary_Datasets.xlsx", sheet = "Supplementary_Data_3") %>%
  rename(Trait=label)

dat1$Supercluster <- factor(dat1$Supercluster, 
                            levels=rev(order_ct))
dat1$Trait <- factor(dat1$Trait, 
                    levels=order_final_traits)
#- bubble plot
p <- ggplot(dat1,
            aes(x=Trait, y=Supercluster,
                fill = -log10(P), 
                size = if.sig.fdr,
                color = if.sig.fdr)) + 
  geom_point(alpha = .85, shape = 21) +
  theme(panel.background = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle = 30, hjust=1),
        legend.key = element_rect(fill = "transparent")) +
  scico::scale_fill_scico(palette = "lajolla",direction=-1, midpoint = 5,
                          name=expression('-log'[10]*'(P)')) + 
  facet_grid(~fct_relevel(disorder_type,
                          "psychiatric",
                          "brain traits",
                          
                          ), 
             scales = "free_x", 
             space = "free_x") +
  scale_size_manual(values=c(1.5, 5))+
  scale_color_manual(values=c("no"="grey","yes"="black")) + 
  labs(size = "Sig.(FDR)")+ ylab("") +xlab("") +
  guides(color="none") + 
  guides(size=guide_legend(override.aes=list(color=c("lightgrey","black")))) 

#- bar plot
dat2 <- read_xlsx("../data/Supplementary_Datasets.xlsx",sheet = "Supplementary_Data_1") 
dat2$label = factor(dat2$label, 
                    levels=order_final_traits)
p1 <- ggplot(dat2 %>% filter(label%in%dat1$Trait),
             aes(x=label, y=ldsc_h2)) +
  geom_bar(stat="identity",width = 1, color="white")+ 
  facet_grid(~fct_relevel(disorder_type,
                          "psychiatric",
                          "brain traits",
                          ), 
             scales = "free_x", 
             space = "free_x") +
  theme(panel.background = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle = 30, hjust=1),
        axis.ticks.x = element_blank(),
        legend.key = element_rect(fill = "transparent")) +
  xlab("")


#- output
ggsave(p,
       filename = "../figures/2a_bubble.pdf",
       width = 16, 
       height = 7.5)
ggsave(p1 + theme(axis.text.x=element_blank()),
       filename = "../figures/2a_bar.pdf",
       width = 16, 
       height = 1.5)

#--- end ---#