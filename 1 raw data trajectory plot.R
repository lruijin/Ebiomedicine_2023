### This is to plot the trajectory paper along all visits i.e., the trajectories
### using z-score
rm(list = ls())
library(reshape2)
library(dplyr)
lipid <- read.csv("../Data/lipid_imputed_wosmoke_combat.csv")
anno <- read.csv("../Data/Lipidomics_annotation_complete.csv",header=T, sep=",")
met_idx = 4:331

lipid_z = lipid
lipid_z[,met_idx] = apply(lipid[,met_idx],2, scale)
lipid_z = lipid_z[,c(1,3,met_idx)]
pca = prcomp(lipid_z[,-c(1,2)])
met.order = data.frame(variable = names(pca$rotation[,1]), 
                       pc1.loading = pca$rotation[,1]) %>%
    mutate(variable = gsub("^.{0,3}", "" , variable)) %>%
    rename(met_ID = variable)
anno = met.order %>% left_join(anno %>%
                                   select(met_ID, Lipid.class))

anno <- anno %>% arrange(Lipid.class, pc1.loading) %>%
    mutate(order = row_number())

lipid_z = melt(lipid_z, id.vars = c("ID", "Visit")) 
names(lipid_z) = c("ID", "Visit", "met_ID", "Lipid")

lipid_z$met_ID = gsub("^.{0,3}", "" , lipid_z$met_ID)
lipid_z = lipid_z %>% 
    left_join(anno %>% 
                  select(met_ID, Lipid.class, pc1.loading, order),
              by = "met_ID")

#### create annotations for the y-axis

# lipid_z$Annotation_comp <- factor(lipid_z$Annotation_comp, 
#     levels = levels(forcats::fct_reorder(lipid_z$Annotation_comp, 
#                                          as.numeric(lipid_z$Lipid.class))))

lipid_z$Lipid.class = factor(lipid_z$Lipid.class)
lipid_z$order = factor(lipid_z$order)
lipid_z$Visit = factor(lipid_z$Visit, levels = c(0,1,2,4),
                       labels = c("10 - 14 wks", "15 - 26 wks", "23 - 31 wks", "33-39 wks"))
breaks <- unique(lipid_z %>% select(order, Lipid.class)) %>%
    group_by(Lipid.class) %>%
    arrange(order, .by_group = T) %>%
    filter(row_number() == ceiling(n()/2))

p = ggplot(lipid_z, aes(x=as.factor(ID), y = order))+
        geom_tile(aes(fill = Lipid))+
        scale_fill_gradient2(name = "Z-scores", 
                             midpoint = 0, low = "blue",mid = "white",high = "red", 
                             limits = c(-3,3), breaks = c(-2,0,2)) + 
        scale_y_discrete(breaks = breaks$order,
                         labels = breaks$Lipid.class)+
        facet_grid(Lipid.class~Visit, switch = "y", scales = "free",
                   space = "free") +
        theme_classic()+
        xlab("Gestational Weeks") +
        ylab("")+
        theme(panel.spacing.y = unit(0.2, "lines"),
              strip.text.y.left = element_blank(),
              axis.text.x.bottom = element_blank(),
              text = element_text(size = 22))
ggsave("../Figures/draft/new figures for trajectories/lipid_long_new.png",
       plot = p, width = 35, height = 42, units = "cm")
