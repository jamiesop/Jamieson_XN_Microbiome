
#
# Title: Figures for XMaSI Trial - targeted metabolites and microbiome 
# Author: Paige Jamieson 
# 

# Environment __________________________________________________________________

library(tidyverse)
library(phyloseq)
library(cowplot)
library(pheatmap)
library(ggpubr) #For legend extraction
library(ggtext) #Applies markdown text elements to ggplot

## Helper Functions

# Standard Error Function
ser <- function(x){
        sd(x)/sqrt(length(x))
}


## Load Data Sets
xns_f <- read.csv("~/Projects/XMaSI/Data/xmas1_xn_fecal_metabs_ng_per_g.csv") #%>% 
#Convert fecal concentration to ?g/g (plasma and urine remain as ng/mL)
#mutate_at(vars(starts_with('f_')), ~.x/1000)

xns_p <- read.csv('~/Projects/XMaSI/Data/XN_metabolites_plasma_ng_per_mL.csv')

xns_u <- read.csv('~/Projects/XMaSI/Data/XN_metabolites_urine_ng_per_mL.csv')

bas <- read.csv('~/Projects/XMaSI/Data/xmas1_bileacids_ug_per_g.csv') %>%
        column_to_rownames('id') %>% 
        rename_with(~paste0('BA_', .x)) %>% 
        mutate(BA_T = dplyr::select(., starts_with("BA_")) %>%  rowSums(na.rm = T),
               BA_CP = dplyr::select(., c('BA_GCA', 'BA_TCA', 'BA_GCDCA', 'BA_TCDCA')) %>% rowSums(na.rm = T),
               BA_GCP = dplyr::select(., c('BA_GCA', 'BA_GCDCA')) %>% rowSums(na.rm = T),
               BA_AGC = dplyr::select(., starts_with('BA_G')) %>% rowSums(na.rm = T),
               BA_UCP = dplyr::select(., c('BA_CA', 'BA_CDCA')) %>% rowSums(na.rm = T),
               BA_AUC = dplyr::select(., -starts_with(c('BA_G', 'BA_T'))) %>% rowSums(na.rm = T),
               #BA_UC = dplyr::select(., c('BA_CA', 'BA_CDCA', 'BA_isoLCA', 'BA_DHLCA', 'BA_NDCA', 'BA_X7.KCDCA', 'BA_UDCA', 'BA_HDCA', 'BA_DHCA', 'BA_X12.KCDCA', 'BA_X7.KDCA', 'BA_UCA', 'BA_MCA')) %>% rowSums(na.rm = T),
               BA_UC = dplyr::select(., c('BA_CA', 'BA_CDCA', 'BA_isoLCA', 'BA_DHLCA', 'BA_NDCA', 'BA_X7.KCDCA', 'BA_UDCA', 'BA_HDCA', 'BA_DHCA', 'BA_X12.KCDCA', 'BA_X7.KDCA', 'BA_UCA', 'BA_MCA')) %>% rowSums(na.rm = T),
               BA_CS = dplyr::select(., c('BA_GDCA', 'BA_TDCA', 'BA_GLCA', 'BA_TLCA')) %>% rowSums(na.rm = T),
               BA_UCS = dplyr::select(., c('BA_LCA', 'BA_DCA')) %>% rowSums(na.rm = T)) %>% 
        rownames_to_column('id')

scfas <- read.csv('~/R/Data/xmas1/xmas1_scfa_mg_per_g.csv') %>% 
        rename('2-methylbutyric' = 'X2.methylbutyric') %>% 
        column_to_rownames('id') %>% 
        rename_with(~paste0('FA_', .x)) %>% 
        rownames_to_column('id')

metas <- read.csv('~/Projects/XMaSI/Data/xmas1_metadata.csv') %>% 
        mutate(time_cont = ifelse(time == "T1", 0,
                                  ifelse(time == "T2", 14, 
                                         ifelse(time == "T3", 28,
                                                ifelse(time == "T4", 42, 56))))) %>% 
        dplyr::rename(BMI = "BMI..kg.m2.") %>% 
        mutate(BMI_cat = ifelse(BMI >= 24.5, "overweight", "healthy")) %>% 
        modify_at('BMI_cat', as.factor) %>% 
        mutate(age_bin = ifelse(age.scn <= 27, 1,
                                ifelse(age.scn <= 33, 2,
                                       ifelse(age.scn <= 38, 3, 4)))) %>% 
        modify_at('age_bin', as.factor) %>% 
        mutate_at(vars(study.id), ~ str_replace(., "XMaS", "")) %>% 
        dplyr::select(-c(specimen.type, date.collected, storage.location, birth.date)) %>% 
        filter(!study.id %in% c('120', '125', '128')) %>% 
        modify_at('study.id', as.factor)

fullxns <- xns_f %>% 
        left_join(., xns_p) %>% 
        left_join(., xns_u) %>% 
        left_join(., metas) %>% 
        filter(treatment == "B") %>% 
        #Convert fecal concentration to ?g/g (plasma and urine remain as ng/mL)
        mutate_at(vars(starts_with('f_')), ~.x/1000)


# load original Phyloseq object
ps_one <- readRDS("~/Projects/XMaSI/Data/xmas1_ps.rds")
# Change meta data in phyloseq object
tax <- tax_table(ps_one)
asv <- otu_table(ps_one)
meta <- metas %>% column_to_rownames('id') %>% 
        sample_data()
gms <- phyloseq(otu_table(asv, taxa_are_rows = FALSE),
                sample_data(meta),
                tax_table(tax))

#Agglomerate to the genus level
ps_genera <- gms %>% tax_glom(taxrank = "Genus") # number of ASVs = 1056
#Remove taxa not seen more than 3 times in at least 20% of the samples
ps_counts <- ps_genera %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE) # number of ASVs = 132
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) )
#Filter out low abundance (>1e-5) taxa >>> does not change number of ASVs
#ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE)


# Main Figures  ----------------------------------------------------------------

# Figure 2. Enterotype Clustering ----------------------------------------------

# 2A [3000w x 3000h]
# 'clusters' variable (dataframe) is from analysis script 

# Plot ordinations with clusters
clus_pal <- c("#f26c29", "#f0cd07", "#5afaba")
names(clus_pal) <- c("1", "2", "3")

ggplot(clusters, aes(x = Axis.1, y = Axis.2, color = clus)) +
        geom_point(size = 18) +
        xlab('Axis 1 [72.5%]') +
        ylab('Axis 2 [15.5%]') +
        scale_y_continuous(limits = c(-0.15, 0.15), breaks = c(-0.15, -0.10, -0.05, 0.00, 0.05, 0.10,0.15)) +
        scale_x_continuous(limits = c(-0.4, 0.3)) +
        theme_cowplot() + 
        theme(axis.title.x = element_text(size=80, vjust = -1),
              axis.title.y = element_text(size=80),
              axis.text.x = element_text(size = 70, vjust = - 1),
              axis.text.y = element_text(size = 70),
              axis.line = element_line(size = 5),
              plot.title = element_text(hjust = 0.5, size = 80),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 60),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=70),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none') +
        scale_color_manual(values = clus_pal, labels = c(expression(''~italic(Prevotella)~' ET'), expression(''~italic(Bacteroides)~' ET'), expression(''~italic(Ruminococcus)~' ET')),
                           name = "Enterotype")

# Extract legend
L <- ggplot(clusters, aes(x = Axis.1, y = Axis.2, color = clus)) +
        geom_point(size = 20) +
        xlab('Axis 1 [72.5%]') +
        ylab('Axis 2 [15.5%]') +
        theme_cowplot() + 
        scale_color_manual(values = clus_pal, labels = c("*Prevotella*-dominant", "*Bacteroides*-dominant", "*Ruminococcus*-dominant"),
                           name = "Enterotype Cluster") +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.title = element_text(size = 55),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.text = element_markdown(size = 55))
legend <- ggpubr::get_legend(L)
ggpubr::as_ggplot(legend)


# 2B [2000w x 6000h]

clus4plot <- clusters %>% 
        dplyr::filter(time == "T1") %>% 
        pivot_longer(cols = c(ASV361, ASV1908, ASV908), names_to = "ET", values_to = "abund")

ggplot(data = clus4plot, aes(x = clus, y = abund, color = clus)) + 
        geom_boxplot(show.legend = F, linewidth = 10) + 
        #geom_point(position = position_jitterdodge(), alpha = 0.6, size = 2) +
        facet_wrap(~ ET, ncol = 1, scales = "free", 
                   labeller = labeller(ET = c("ASV1908" = "Bacteroides", "ASV361" = "Prevotella", "ASV908" = "Ruminococcus"))) +
        cowplot::theme_cowplot() + 
        scale_color_manual(values = clus_pal, labels = c("*Prevotella*", "*Bacteroides*", "*Ruminococcus*"),
                           name = "Enterotype") +
        #scale_x_discrete(labels = c('*Prev*-dominant', '*Bact*-dominant', '*Rumin*-dominant'),
                         #guide = guide_axis(angle = 45)) +
        #xlab("Cluster") +
        ylab("Taxa Abundance at Baseline") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size=120),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 10),
              plot.title = element_text(hjust = 0.5, size = 90, face = 'italic'),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=95, face = 'italic'),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

# 2C
# NB model of differential taxa 
clus.nb <- clusters %>% 
        #simplify model with 3 time points instead of 5
        dplyr::filter(time %in% c('T1', 'T3', 'T5')) %>% 
        # Remove ASV's with too many zeros, not enough information for model
        dplyr::select(-c(ASV14, ASV647, ASV658, ASV1146, ASV1520, ASV1847, ASV1949, ASV2121)) %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'value') %>% 
        group_by(ASV) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) glmmTMB(value ~ clus*time + (1|study.id), family=nbinom2, data = x))) %>% 
        # test contrast for each cluster comparing taxa at close-out to taxa at baseline
        mutate(clus1_t5vt1 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,0,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(clus2_t5vt1 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,1,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(clus3_t5vt1 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,0,1), nrow = 1)))[[1,6]][1])) %>% 
        ungroup() %>% 
        mutate(clus1_adj = p.adjust(clus1_t5vt1, method = 'BH')) %>% 
        mutate(clus2_adj = p.adjust(clus2_t5vt1, method = 'BH')) %>%
        mutate(clus3_adj = p.adjust(clus3_t5vt1, method = 'BH'))

# filter for significant differentially abudnant taxa 
clus.filtered <- clus.nb %>% 
        dplyr::filter(clus1_adj <= 0.05 | clus2_adj <= 0.05) %>% 
        left_join(., alltax) %>% 
        pivot_longer(cols = ends_with("_adj"), names_to = "cluster", values_to = "padj") %>% 
        mutate(clus = ifelse(cluster == "clus1_adj", "1", 
                             ifelse(cluster == "clus2_adj", 2, 3))) %>% 
        mutate(clus_ASV = paste0(clus, "_", ASV))

# Log2 Fold Change Function
l2fc <- function(x, condB, condA){
        cB <- x[[condB]] %>% mean()
        cA <- x[[condA]] %>% mean()
        l2fc <- log2(cB/cA)
        return(l2fc)
}

# Pull differentially abundant ASV's
db_asv <- clus.filtered$ASV

#Filter count data by asv's and T1, T5, pivot long time points and ASVs, 
fcData <- clusters %>% 
        dplyr::filter(time %in% c("T1", "T5")) %>% 
        dplyr::select(all_of(db_asv), study.id, clus, time) %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'count') %>% 
        group_by(clus, ASV) %>% 
        nest() %>% 
        mutate(data = purrr::map(data, function(x) pivot_wider(x, names_from = "time", values_from = "count"))) %>% 
        mutate(l2fc = purrr::map_dbl(data, function(x) l2fc(x, condB="T5", condA="T1"))) %>% 
        mutate_at(vars(l2fc), ~replace(., is.nan(.), 0)) %>% 
        left_join(., alltax, by = "ASV") %>% 
        mutate(clus_ASV = paste0(clus, "_", ASV)) %>% 
        left_join(clus.filtered %>% dplyr::select(clus_ASV, padj))

#Define pallette colors
algae = colorRampPalette(c("#1A286E", "#E5FA7B"))(100)

# Create dataframe for heatmap values
m1 <- fcData %>% 
        dplyr::select(clus, l2fc, ASV) %>% 
        pivot_wider(names_from = "clus", values_from = "l2fc") %>% 
        column_to_rownames('ASV') 
# Create dataframe for heatmap significance
m2 <- fcData %>% 
        dplyr::select(clus, padj, ASV) %>% 
        mutate(sig = ifelse(padj <= 0.001, "***",
                            ifelse(padj <= 0.01, "**",
                                   ifelse(padj <= 0.05, "*", " ")))) %>% 
        dplyr::select(-padj) %>% 
        pivot_wider(names_from = "clus", values_from = c("sig")) %>% 
        column_to_rownames('ASV') 


#Define labels for 'ASV' and clusters
asv_labs <- c("Blautia", "f_Lachnospiraceae_ASV64",            
              "Roseburia", "f_Lachnospiraceae_ASV501",            
              "f_Oscillospiraceae_ASV683", "Flavonifractor",                     
              "Dialister", "o_Clostridia_vadinBB60_group_ASV1662")

clus_labs <- c("Prevotella (1)", "Bacteroides (2)", "Ruminococcus (3)")        


pheatmap::pheatmap(m1, color = algae, cluster_cols = F, cluster_rows = F,
                   display_numbers = m2, fontsize_number=20, cellheight=20, number_color = "white",
                   labels_cols = asv_labs)


# Heatmap
pheatmap::pheatmap(m1, color = algae, cluster_cols = F, cluster_rows = F,
                   display_numbers = m2, fontsize_number=18, cellheight=25, number_color = "white",
                   labels_row = asv_labs, angle_col = 0,
                   labels_col = clus_labs,
                   fontsize_row = 12,
                   fontsize_col = 12, 
                   cellwidth = 25)


# Figure 3. XN Metabolism ------------------------------------------------------

metab_pal <- c('#F8766D', '#CF6339', '#D89000', '#F2C634', '#A3A500', '#ADDC3D', '#39B600', '#00BF7D', '#00BFC4', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC', '#BF2256')

names(metab_pal) <- c('101', '102', '103','105','111','113','115','116','123','124','127','129','130','131')

xns_summary <- fullxns %>% 
        pivot_longer(cols = c(XN_f_XN:XN_u_DDXN), names_to = 'metab', values_to = 'conc') %>% 
        group_by(time_cont, metab) %>% 
        summarise(mean = mean(conc)) %>% 
        pivot_wider(names_from = 'metab', values_from = 'mean')

#Plot of XN concentrations
ggplot(fullxns, aes(x=time_cont, y = f_xn, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(linewidth = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'XN (?g/g dry weight)', limits = c(0, 180), breaks = seq(0,180, by = 20)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(linewidth = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of IXN concentrations
ggplot(fullxns, aes(x=time_cont, y = f_ixn, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'IXN (?g/g dry weight)', limits = c(0, 90), breaks = seq(0,90, by = 10)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of 8PN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_f_X8.PN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '8PN (?g/g dry weight)', limits = c(0, 12), breaks = seq(0,12, by = 2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of DXN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_f_DXN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'DXN (?g/g dry weight)', limits = c(0, 8), breaks = seq(0,8, by = 1)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of 6PN concentrations
ggplot(fullxns, aes(x=time_cont, y = f_6pn, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '6PN (?g/g dry weight)', limits = c(0, 1), breaks = seq(0,1, by = 0.2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of DDXN concentrations
ggplot(fullxns, aes(x=time_cont, y = f_ddxn, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(linewidth = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'DDXN (?g/g dry weight)', limits = c(0, 4), breaks = seq(0,4, by = 1)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of total metabolite concentrations
ggplot(fullxns, aes(x=time_cont, y = f_total, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(linewidth = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'total XN metabolites (?g/g dry weight)', limits = c(0, 180), breaks = seq(0,180, by = 20)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(linewidth = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')


#Extract legend (save as 500w x 2000h)

p1 <- ggplot(fullxns, aes(x=time_cont, y = XN_f_X6.PN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 30) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '6PN (?g/g dry weight)', limits = c(0, 1.2), breaks = seq(0,1.2, by = 0.2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"))

legend <- ggpubr::get_legend(p1)

as_ggplot(legend)


# Figure 4. & 5. DXN & DDXN Spectra (no code) ----------------------------------


# Figure 6. XN Metabolism by Enterotype [1800w x 2200h] ------------------------

clus1 <- clusters %>% 
        # change metabolite concentration from ng/g to ug/g
        dplyr::mutate_at(vars(starts_with('f_')), ~.x/1000) %>% 
        dplyr::select(-c(f_xn, f_ixn)) %>% 
        dplyr::filter(base_clus == '1' & time != 'T1') %>% 
        pivot_longer(starts_with('f_'), names_to = 'metab', values_to = 'value')
clus2 <- clusters %>% 
        # change metabolite concentration from ng/g to ug/g
        dplyr::mutate_at(vars(starts_with('f_')), ~.x/1000) %>% 
        dplyr::select(-c(f_xn, f_ixn)) %>%
        dplyr::filter(base_clus == '2' & time != 'T1') %>% 
        pivot_longer(starts_with('f_'), names_to = 'metab', values_to = 'value')
clus3 <- clusters %>% 
        # change metabolite concentration from ng/g to ug/g
        dplyr::mutate_at(vars(starts_with('f_')), ~.x/1000) %>% 
        dplyr::select(-c(f_xn, f_ixn)) %>%
        dplyr::filter(base_clus == '3' & time != 'T1') %>% 
        pivot_longer(starts_with('f_'), names_to = 'metab', values_to = 'value')
        
library(RColorBrewer)

colors4metabs = c('#EF894D', '#F35466', '#7FCF6F', '#66A0C9')
names(colors4metabs) = c('6PN', '8PN', 'DDXN', 'DXN')


colourCount = length(unique(clus1$metab))
getPalette = colorRampPalette(brewer.pal(14, "Paired"))

ggplot(clus1, aes(x = time, y = value, fill = metab)) +
        geom_bar(stat = 'identity', size = 0, show.legend = F) +
        scale_x_discrete(labels = c('14', '28', '42', '56'))+
        scale_y_continuous(limits = c(0,35)) +
        scale_fill_manual(values=colors4metabs,
                          name = 'Metabolites', labels = c('6PN', '8PN', 'DDXN', 'DXN', 'IXN'),
                          ) + 
        guides(col = guide_legend(ncol = 1)) +
        theme(axis.text.x = element_text(angle = 0)) +
        cowplot::theme_cowplot() +
        xlab("Time (Days)") +
        ylab("Concentration (µg/g dry wt)") +
        ggtitle(expression(''~italic(Prevotella)~' ET')) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(linewidth = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

ggplot(clus2, aes(x = time, y = value, fill = metab)) +
        geom_bar(stat = 'identity', size = 0, show.legend = F) +
        scale_x_discrete(labels = c('14', '28', '42', '56'))+
        scale_y_continuous(limits = c(0,35)) +
        scale_fill_manual(values=colors4metabs,
                          name = 'Metabolites', labels = c('6PN', '8PN', 'DDXN', 'DXN', 'IXN')) +  
        guides(col = guide_legend(ncol = 1)) +
        theme(axis.text.x = element_text(angle = 0)) +
        cowplot::theme_cowplot() +
        xlab("Time (Days)") +
        ylab("Concentration (µg/g dry wt)") +
        ggtitle(expression(''~italic(Bacteroides)~' ET')) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(linewidth = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

ggplot(clus3, aes(x = time, y = value, fill = metab)) +
        geom_bar(stat = 'identity', size = 0, show.legend = F) +
        scale_x_discrete(labels = c('14', '28', '42', '56'))+
        scale_y_continuous(limits = c(0,35)) +
        scale_fill_manual(values=colors4metabs,
                          name = 'Metabolites', labels = c('6PN', '8PN', 'DDXN', 'DXN', 'IXN')) + 
        guides(col = guide_legend(ncol = 1)) +
        theme(axis.text.x = element_text(angle = 0)) +
        cowplot::theme_cowplot() +
        xlab("Time (Days)") +
        ylab("Concentration (µg/g dry wt)") +
        ggtitle(expression(''~italic(Ruminococcus)~' ET')) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(linewidth = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

# Extract legend
L1 <- ggplot(clus3, aes(x = time, y = value, fill = metab)) +
        geom_bar(stat = 'identity', size = 0) +
        scale_x_discrete(labels = c('14', '28', '42', '56'))+
        scale_fill_manual(values=colors4metabs,
                          name = 'Metabolites', labels = c('6PN', '8PN', 'DDXN', 'DXN')) + 
        guides(col = guide_legend(ncol = 1)) +
        theme(axis.text.x = element_text(angle = 0)) +
        cowplot::theme_cowplot() +
        xlab("Time (Days)") +
        ylab("Concentration (µg/g)") +
        ggtitle(expression(''~italic(Ruminococcus)~' ET')) +
        theme(
              legend.text = element_text(size = 60),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=40),
              legend.key.size = unit(2, "cm"))
legend <- ggpubr::get_legend(L1)
ggpubr::as_ggplot(legend)


# Figure 7. 8PN by BMI ---------------------------------------------------------

xn_by_bmi_fig <- fullxns %>% 
        filter(treatment == 'B') %>% 
        pivot_longer(cols = c('XN_f_X8.PN','XN_p_X8.PN', 'XN_u_X8.PN'),
                     names_to = 'metabolite') %>% 
        group_by(BMI_cat, time_cont, metabolite) %>%
        summarise(mean = mean(value),
                  SEM = ser(value))

#Set pallete
bmi_pal <- c('#368996', '#B82E00')
bmi_pal <- c('#18C591', '#F6A94B')
names(bmi_pal) <- c('healthy', 'overweight')
#shapes
bmi_sh <- c(15, 17)
names(bmi_sh) <- c('healthy', 'overweight')

#Subset by sample site
x8pn_ur_data <- xn_by_bmi_fig %>% 
        filter(metabolite == 'XN_u_X8.PN')
x8pn_pl_data <- xn_by_bmi_fig %>% 
        filter(metabolite == 'XN_p_X8.PN')
x8pn_st_data <- xn_by_bmi_fig %>% 
        filter(metabolite == 'XN_f_X8.PN')

#Fecal 8PN by BMI
ggplot(x8pn_st_data, aes(x = time_cont, y = mean, color = BMI_cat, shape = BMI_cat))  +
        geom_point(size = 12) +
        geom_path(aes(group = BMI_cat), size = 5) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), width = 2, size = 5) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = bmi_pal, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_shape_manual(values = bmi_sh, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = 'Fecal 8PN', x = 'Time (Days)', y = 'Fecal 8PN (?g/g dry weight)') +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              legend.position = 'none')

#Plasma 8PN by BMI
ggplot(x8pn_pl_data, aes(x = time_cont, y = mean, color = BMI_cat, shape = BMI_cat))  +
        geom_point(size = 12) +
        geom_path(aes(group = BMI_cat), size = 5) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), width = 2, size = 5) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = bmi_pal, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_shape_manual(values = bmi_sh, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = 'Plasma 8PN', x = 'Time (Days)', y = 'Plasma 8PN (ng/mL)') +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              legend.position = 'none')

#Urine 8PN by BMI
ggplot(x8pn_ur_data, aes(x = time_cont, y = mean, color = BMI_cat, shape = BMI_cat))  +
        geom_point(size = 12) +
        geom_path(aes(group = BMI_cat), size = 5) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), width = 2, size = 5) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = bmi_pal, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_shape_manual(values = bmi_sh, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = '24-hr Urine 8PN', x = 'Time (Days)', y = '24-hr Urine 8PN (ng/mL)') +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              legend.position = 'none')

#Extract legend
p1 <- ggplot(x8pn_ur_data, aes(x = time_cont, y = mean, color = BMI_cat, shape = BMI_cat))  +
        geom_point(size = 20) +
        geom_path(show.legend = FALSE, size = 5) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = bmi_pal, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_shape_manual(values = bmi_sh, labels = c('normal', 'overweight'),
                           name = 'BMI') +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = '24-hr Urine 8PN', x = 'Time (Days)', y = '24-hr Urine 8PN (ng/mL)') +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80))
legend <- ggpubr::get_legend(p1)
as_ggplot(legend)


# Figure 8. sPLS Heatmap of Microbe by Metabolite Associations (pulled from MixOmics analysis)


# Figure 9. Bile Acids [2000w x 1800h] -----------------------------------------

# Summarize Bile Acid groups
ba_groups4plot <- metas %>% 
        left_join(bas) %>% 
        dplyr::select(id, treatment, time, time_cont, study.id, BA_AUC, BA_UCS) %>% 
        # Change ug/g to mg/g BA Groups
        mutate(across(.cols = starts_with('BA_'), ~.x/1000)) %>% 
        pivot_longer(starts_with('BA_'), 
                     names_to = 'metab') %>% 
        group_by(treatment, time_cont, metab) %>% 
        summarise(mean = mean(value),
                  SEM = ser(value))

# Extract Bile Acid groups found significantly different 
AUCgroup_data <- ba_groups4plot %>% 
        dplyr::filter(metab == 'BA_AUC')
UCSgroup_data <- ba_groups4plot %>% 
        dplyr::filter(metab == 'BA_UCS')

ba_pal <- c('#EF115F', '#1AD0BE')
names(ba_pal) <- c('A', 'B')

#Unconjugated BA 
ggplot(data = AUCgroup_data, aes(x = time_cont, y = mean, group = treatment, color = treatment))+
        geom_point(size = 12, position = pd, show.legend = F) +
        geom_path(size = 5, show.legend = F, position = pd) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = pd, width = 5, size = 5, show.legend = F) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = ba_pal, labels = c('Placebo', 'XN-treated'),
                           name = 'Group') +
        scale_y_continuous(limits = c(10, 35)) +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = "Unconjugated Bile Acids", x = 'Time (Days)', y = 'Concentration (mg/g dry wt)') +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 60),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Secondary BA 
ggplot(data = UCSgroup_data, aes(x = time_cont, y = mean, group = treatment, color = treatment))+
        geom_point(size = 12, position = pd, show.legend = F) +
        geom_path(size = 5, show.legend = F, position = pd) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = pd, width = 5, size = 5, show.legend = F) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = ba_pal, labels = c('Placebo', 'XN-treated'),
                           name = 'Group') +
        scale_y_continuous(limits = c(6, 18)) +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = "Secondary Bile Acids", x = 'Time (Days)', y = 'Concentration (mg/g dry wt)') +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 60),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

# Extract Legend [Use default size]
L <- ggplot(data = UCSgroup_data, aes(x = time_cont, y = mean, group = treatment, color = treatment))+
        geom_point(size = 20, position = pd) +
        geom_path(size = 5, show.legend = F, position = pd) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = pd, width = 5, size = 5, show.legend = F) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = ba_pal, labels = c('Placebo', 'XN-treated'),
                           name = 'Group') +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = "Unconjugated Bile Acids", x = 'Time (Days)', y = 'Concentration (mg/g dry weight)') +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.title = element_text(size = 55),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.text = element_markdown(size = 55))
legend <- ggpubr::get_legend(L)
ggpubr::as_ggplot(legend)


# Figure 10. Bile Acids by Enterotype [2000w x 1800h] ---------------------------

ba4plot <- clusters %>% 
        rownames_to_column('id') %>% 
        left_join(bas) %>% 
        dplyr::select(id, base_clus, time, time_cont, study.id, BA_AUC, BA_UCS) %>% 
        # Change ug/g to mg/g BA Groups
        mutate(across(.cols = starts_with('BA_'), ~.x/1000)) %>% 
        pivot_longer(starts_with('BA_'), 
                     names_to = 'metab') %>% 
        group_by(base_clus, time_cont, metab) %>% 
        summarise(mean = mean(value),
                  SEM = ser(value))

AUC_data <- ba4plot %>% 
        dplyr::filter(metab == 'BA_AUC')
UCS_data <- ba4plot %>% 
        dplyr::filter(metab == 'BA_UCS')

#cluster palette
clus_pal <- c("#f26c29", "#f0cd07", "#5afaba")
names(clus_pal) <- c("1", "2", "3")

# to offset errorbars
pd <- position_dodge(width = 2.5)

ggplot(AUC_data, aes(x = time_cont, y = mean, color = base_clus))  +
        geom_point(size = 12, position = pd, show.legend = F) +
        geom_path(aes(group = base_clus), size = 5, show.legend = F) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = pd, width = 5, size = 5, show.legend = F) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = clus_pal, labels = c("Prevotella", "Bacteroides", "Ruminococcus"),
                           name = "Enterotype") +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = "Unconjugated Bile Acids", x = 'Time (Days)', y = 'Concentration (mg/g dry wt)') +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 60),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')
        
ggplot(UCS_data, aes(x = time_cont, y = mean, color = base_clus))  +
        geom_point(size = 12, position = pd, show.legend = F) +
        geom_path(aes(group = base_clus), size = 5, show.legend = F) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = pd, width = 5, size = 5, show.legend = F) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = clus_pal, labels = c("Prevotella", "Bacteroides", "Ruminococcus"),
                           name = "Enterotype") +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = "Secondary Bile Acids", x = 'Time (Days)', y = 'Concentration (mg/g dry wt)') +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 60),
              legend.title = element_text(size = 60),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

# Extract Legend [Use Default Size]
L <- ggplot(UCS_data, aes(x = time_cont, y = mean, color = base_clus))  +
        geom_point(size = 20, position = pd) +
        geom_path(aes(group = base_clus), size = 5, show.legend = F) +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = pd, width = 5, size = 5, show.legend = F) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = clus_pal, labels = c("*Prevotella*", "*Bacteroides*", "*Ruminococcus*"),
                           name = "Enterotype") +
        scale_x_continuous(breaks = c(0,14,28,42,56), labels = c(0,14,28,42,56)) +
        labs(title = "Secondary Bile Acids", x = 'Time (Days)', y = 'Concentration (mg/g dry wt)') +
        theme(axis.title.x = element_text(size=70, vjust = -1),
              axis.title.y = element_text(size=70),
              axis.text.x = element_text(size = 60, vjust = - 1),
              axis.text.y = element_text(size = 60),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 70),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.title = element_text(size = 55),
              strip.text = element_text(size=60),
              panel.spacing = unit(5, "lines"),
              legend.text = element_markdown(size = 55))
legend <- ggpubr::get_legend(L)
ggpubr::as_ggplot(legend)


# Supplemental Figures  --------------------------------------------------------

# Figure S1. CONSORT Diagram (no code) -----------------------------------------


# Figure S2. Alpha-diversity [5000w x 2000h] -----------------------------------

# Calculate alpha diversity measures
adiv <- gms %>% 
        estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
        rownames_to_column("id") %>% 
        mutate_at(vars(id), ~ str_replace(., fixed("."), "-"))

adiv_pal <- c('#368996', '#B82E00')
names(adiv_pal) <- c('B', 'A')

adiv4plot <- adiv %>% left_join(., data.frame(sample_data(gms)) %>% rownames_to_column('id')) %>% 
        pivot_longer(., cols = c(Observed, Shannon, Simpson), names_to = 'measure') 

# Boxplot
ggplot(data = adiv4plot, aes(x = factor(time_cont), y = value, color = treatment)) +
        geom_boxplot(show.legend = F, linewidth = 5) +
        geom_point(position = position_jitterdodge(), alpha = 0.6, size = 12) +
        #geom_jitter(position = position_dodge(0.6-0.8)) +
        facet_wrap(~measure, ncol = 3, scales = 'free') +
        cowplot::theme_cowplot()+
        scale_color_manual(values = adiv_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        xlab('Time (Days)') +
        ylab('')+
        scale_x_discrete(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"))


# Figure S4. XN Metabolism (plasma) [2000w x 2000h] ----------------------------

#Plot of XN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_p_XN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'XN (ng/mL)', limits = c(0, 18), breaks = seq(0,18, by = 2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of IXN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_p_IXN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'IXN (ng/mL)', limits = c(0, 18), breaks = seq(0,18, by = 2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of 8PN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_p_X8.PN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '8PN (ng/mL)', limits = c(0, 18), breaks = seq(0,18, by = 2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of DXN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_p_DXN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'DXN (ng/mL)', limits = c(0, 18), breaks = seq(0,18, by = 2)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of 6PN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_p_X6.PN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '6PN (ng/mL)', limits = c(0, 18), breaks = seq(0,18, by = 12)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')


# Figure S5. XN Metabolism (urine) [2000w x 2000h] -----------------------------

#Plot of XN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_u_XN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'XN (ng/mL)', limits = c(0, 60), breaks = seq(0,60, by = 5)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of IXN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_u_IXN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'IXN (ng/mL)', limits = c(0, 300), breaks = seq(0,300, by = 50)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of 8PN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_u_X8.PN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '8PN (ng/mL)', limits = c(0, 60), breaks = seq(0,60, by = 5)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of DXN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_u_DXN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'DXN (ng/mL)', limits = c(0, 25), breaks = seq(0,25, by = 5)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of 6PN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_u_X6.PN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= '6PN (ng/mL)', limits = c(0, 25), breaks = seq(0,25, by = 5)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')

#Plot of DDXN concentrations
ggplot(fullxns, aes(x=time_cont, y = XN_u_DDXN, color = study.id, group = study.id)) +
        geom_bar(aes(group = time_cont), stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
        geom_path(size = 4.5, show.legend = FALSE) +
        geom_point(size = 9) + 
        scale_color_manual(name = 'Subject ID', 
                           values = metab_pal) +
        theme_minimal_hgrid() +
        scale_x_continuous(name='Time (Days)', breaks=c(0,14,28,42,56), expand = c(0.005, 0.005)) +
        scale_y_continuous(name= 'DDXN (ng/mL)', limits = c(0, 5), breaks = seq(0,5, by = 1)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = 'none')


# Figure S6. Spearman correlation [2000w x 1800h] ------------------------------

# f vs p XN
ggplot(fullxns, aes(x=f_xn, y=p_xn))+
        geom_point(size=8) +
        labs(title = 'XN', x = 'Fecal XN (?g/g dry weight)', y ='Plasma XN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))

# f vs p IXN
ggplot(fullxns, aes(x=f_ixn, y=p_ixn))+
        geom_point(size=8)+
        labs(title = 'IXN', x = 'Fecal IXN (?g/g dry weight)', y ='Plasma IXN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm')) 

# f vs p 8PN
ggplot(fullxns, aes(x=f_8pn, y=p_8pn))+
        geom_point(size=8)+
        labs(title = '8PN', x = 'Fecal 8PN (?g/g dry weight)', y ='Plasma 8PN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))

# f vs p DXN
ggplot(fullxns, aes(x=f_dxn, y=p_dxn))+
        geom_point(size=8)+
        labs(title = 'DXN', x = 'Fecal DXN (?g/g dry weight)', y ='Plasma DXN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm')) 

# f vs p DDXN
ggplot(fullxns, aes(x=f_ddxn, y=XN_p_DDXN))+
        geom_point(size=8)+
        labs(title = 'DDXN', x = 'Fecal DDXN (?g/g dry weight)', y ='Plasma DDXN (ng/mL)')+
        theme_cowplot(12) + 
        scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by = 0.25)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))


# f vs u XN
ggplot(fullxns, aes(x=XN_f_XN, y=XN_u_XN))+
        geom_point(size=8)+
        labs(title = 'XN', x = 'Fecal XN (?g/g dry weight)', y ='24-hr Urine XN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))

# f vs u IXN
ggplot(fullxns, aes(x=XN_f_IXN, y=XN_u_IXN))+
        geom_point(size=8)+
        labs(title = 'IXN', x = 'Fecal IXN (?g/g dry weight)', y ='24-hr Urine IXN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))

# f vs u 8PN
ggplot(fullxns, aes(x=XN_f_X8.PN, y=XN_u_X8.PN))+
        geom_point(size=8)+
        labs(title = '8PN', x = 'Fecal 8PN (?g/g dry weight)', y ='24-hr Urine 8PN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))

# f vs u DXN
ggplot(fullxns, aes(x=XN_f_DXN, y=XN_u_DXN))+
        geom_point(size=8)+
        labs(title = 'DXN', x = 'Fecal DXN (?g/g dry weight)', y ='24-hr Urine DXN (ng/mL)')+
        theme_cowplot(12) + 
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))

# f vs u DDXN
ggplot(fullxns, aes(x=f_ddxn, y=XN_u_DDXN))+
        geom_point(size=8)+
        labs(title = 'DDXN', x = 'Fecal DDXN (?g/g dry weight)', y ='24-hr urine DDXN (ng/mL)')+
        theme_cowplot(12) + 
        scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by = 0.25)) +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'))


# Figure S7. Short-chain Fatty Acids -------------------------------------------

#Summarize mean and standard error
fullscfa <- left_join(scfas, metas) %>% 
        pivot_longer(., cols = c(acetic:isovaleric), names_to = 'scfa') %>% 
        group_by(time_cont, scfa, treatment) %>% 
        summarise(mean = mean(value),
                  SEM = ser(value))

#Set color palette
scfa_pal <- c('#368996', '#B82E00')
names(scfa_pal) <- c('B', 'A')

ggplot(data = fullscfa, aes(x = time_cont, y = mean, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        facet_wrap(~scfa, ncol=3, scales = 'free') +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = scfa_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        xlab('Time (Days)') +
        ylab('Concentration (mg/g dry weight)') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        #labs(title = 'Fecal 8PN')
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"))


# Figure S8. Plasma Cytokines --------------------------------------------------

#Summarize mean and standard error
cyto_summ <- left_join(inflam, metas) %>% 
        filter(!study.id %in% c('117', '126')) %>%
        mutate(across(.cols = c(IL12:plasma_LBP), as.numeric)) %>% 
        pivot_longer(., cols = c('IL12', 'IL10', 'TNF'), names_to = 'cytokines') %>% 
        group_by(time_cont, cytokines, treatment) %>% 
        summarise(mean = mean(value),
                  SEM = ser(value))

#Set color palette
cyto_pal <- c('#368996', '#B82E00')
names(cyto_pal) <- c('B', 'A')

ggplot(data = cyto_summ, aes(x = time_cont, y = mean, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        facet_wrap(~cytokines, ncol=3, scales = 'free') +
        geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = cyto_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        xlab('Time (Days)') +
        ylab('Concentration (pg/mL)') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        #labs(title = 'Fecal 8PN')
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"))


# Figure S9. Markers of Gut Barrier --------------------------------------------

gh_summ <- left_join(inflam, metas) %>% 
        filter(!study.id %in% c('117', '126')) %>%
        rename("IFABP" = "plasma_I.FABP") %>% 
        mutate(across(.cols = c(IL12:plasma_LBP), as.numeric)) %>% 
        pivot_longer(., cols = c(plasma_CD14:plasma_LBP), names_to = 'markers') %>% 
        group_by(time_cont, markers, treatment) %>% 
        summarise(mean = mean(value),
                  SEM = ser(value)) %>% 
        pivot_wider(names_from = "markers", values_from = c("mean", "SEM"))

#Set color palette
gh_pal <- c('#368996', '#B82E00')
names(gh_pal) <- c('B', 'A')

#Endotoxin plot (save as 3000w x 3000h)
ggplot(data = gh_summ, aes(x = time_cont, y = mean_endotoxin, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        geom_errorbar(aes(ymin = mean_endotoxin - SEM_endotoxin, ymax = mean_endotoxin + SEM_endotoxin), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = gh_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        labs(title = 'Plasma Endotoxin', x = "Time (Days)", y = "Concentration (EU/mL)") +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = "none")

#Plasma CD14 plot (save as 3000w x 3000h)
ggplot(data = gh_summ, aes(x = time_cont, y = mean_plasma_CD14, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        geom_errorbar(aes(ymin = mean_plasma_CD14 - SEM_plasma_CD14, ymax = mean_plasma_CD14 + SEM_plasma_CD14), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = gh_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        labs(title = 'Plasma soluble CD14', x = "Time (Days)", y = "Concentration (?g/mL)") +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = "none")

#Plasma I-FABP plot (save as 3000w x 3000h)
ggplot(data = gh_summ, aes(x = time_cont, y = mean_IFABP, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        geom_errorbar(aes(ymin = mean_IFABP - SEM_IFABP, ymax = mean_IFABP + SEM_IFABP), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = gh_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        labs(title = 'Plasma I-FABP', x = "Time (Days)", y = "Concentration (pg/mL)") +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = "none")

#Plasma LBP plot (save as 3000w x 3000h)
ggplot(data = gh_summ, aes(x = time_cont, y = mean_plasma_LBP, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        geom_errorbar(aes(ymin = mean_plasma_LBP - SEM_plasma_LBP, ymax = mean_plasma_LBP + SEM_plasma_LBP), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = gh_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        labs(title = 'Plasma LBP', x = "Time (Days)", y = "Concentration (?g/mL)") +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"),
              legend.position = "none")

#Extract legend
p1 <- ggplot(data = gh_summ, aes(x = time_cont, y = mean_plasma_LBP, group = treatment, color = treatment))+
        geom_point(size = 20)+
        geom_path(show.legend = FALSE, size = 5) +
        geom_errorbar(aes(ymin = mean_plasma_LBP - SEM_plasma_LBP, ymax = mean_plasma_LBP + SEM_plasma_LBP), width = 2, size = 5, show.legend = FALSE) +
        cowplot::theme_cowplot() +
        scale_color_manual(values = gh_pal, labels = c('CTR', 'XN'),
                           name = 'Group') +
        scale_x_continuous(breaks = c(0, 14, 28, 42, 56), labels = c(0, 14, 28, 42, 56)) +
        labs(title = 'Plasma LBP', x = "Time (Days)", y = "Concentration (?g/mL)") +
        theme(axis.title.x = element_text(size=90, vjust = -1),
              axis.title.y = element_text(size=90),
              axis.text.x = element_text(size = 80, vjust = - 1),
              axis.text.y = element_text(size = 80),
              axis.line = element_line(size = 6),
              plot.title = element_text(hjust = 0.5, size = 90),
              plot.margin = margin(2,2,2,2, 'cm'),
              legend.text = element_text(size = 80),
              legend.title = element_text(size = 80),
              strip.text = element_text(size=80),
              panel.spacing = unit(5, "lines"))

legend <- ggpubr::get_legend(p1)

as_ggplot(legend)

