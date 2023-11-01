# 
# Gut Microbiota & Microbiota-Metabolism by Enterotypes in the XMaS Trial
# Author: Paige Jamieson
#

# Environment ------------------------------------------------------------------

library(tidyverse) 
library(phyloseq)
library(lmerTest)
library(glmmTMB)
library(multcomp)
library(mixOmics)
library(cluster)

# Helper Functions

# Standard Error
ser <- function(x){
        sd(x)/sqrt(length(x))
}

#Auto Scaling
AS_helper <- function(x) {
        (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(mtb){
        mtb_scaled <- apply(mtb, 2, AS_helper) 
        return(mtb_scaled)
}

# Robust scaling (value - median/IQR)
RS_helper <- function(x){
        (x- median(x)) /(quantile(x,probs = .75)-quantile(x,probs = .25))
}

robust_scale <- function(mtb){
        mtb_scaled <- apply(mtb, 2, RS_helper)
        return(mtb_scaled)
}


# Natural log transform
log_helper <- function(x, min.val){
        log((x + sqrt(x^2 + min.val^2))/2)
}


log_transform <- function(mtb){
        mtb_nz <- mtb[ ,which(apply(mtb, 2, sum) != 0)]
        min.val <- min(abs(mtb_nz[mtb_nz!=0]))/10
        mtb_log_trans <- apply(mtb_nz, 2, log_helper, min.val)
        return(as.data.frame(mtb_log_trans))
}


# Pareto scale
PS_helper <- function(x){
        (x - mean(x))/sqrt(sd(x, na.rm = T))
}

pareto_scale <- function(mtb){
        mtb_scaled <- apply(mtb, 2, PS_helper) 
        return(mtb_scaled)
}


# Center log-ratio functions
geoMean <- function(x) {
        exp(mean(log(x)))
}

centerHelper <- function(x) {
        log(x/geoMean(x))
}

clr <- function(x){
        nz <- data.frame(apply(x, 2, function(y) replace(y, which(y == 0), 1)))
        data.frame(apply(nz, 2, centerHelper))
}

# Load in Data -----------------------------------------------------------------

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

gmarks <- read.csv('~/Projects/XMaSI/Data/xmas1_total_infl_markers.csv') %>% 
        column_to_rownames('id') %>% 
        rename_with(~paste0('GM_', .x)) %>% 
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

# load original Phyloseq object
gms <- readRDS('~/Projects/XMaSI/Data/xmas1_ps.rds')
# number of ASVs = 2152

#Give arbitrary names to the taxa as opposed to keeping as just DNA-sequences which identify them
taxa_names(gms) <- paste0("ASV", seq(ntaxa(gms)))

#Fill in missing genus names with family:
renames <- rownames(tax_table(gms)[is.na(tax_table(gms)[, 'Genus'])])
taxdf <- tax_table(gms)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_table(gms)[renames, 'Genus'] <- renamed_genus

#Fill in missing family names with order:
f_renames <- rownames(tax_table(gms)[is.na(tax_table(gms)[, 'Family'])])
f_taxdf <- tax_table(gms)[f_renames, ]
renamed_family <- unname(sapply(taxa_names(f_taxdf), function(x) paste0('o_', f_taxdf[x, 'Order'], '_', x)))
tax_table(gms)[f_renames, c('Family', 'Genus')] <- renamed_family

#Agglomerate to the genus level
ps_genera <- gms %>% tax_glom(taxrank = "Genus") # number of ASVs = 1056
#Remove taxa not seen more than 3 times in at least 20% of the samples
ps_counts <- ps_genera %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE) # number of ASVs = 132
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) )
#Filter out low abundance (>1e-5) taxa >>> does not change number of ASVs
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE)


# Gut Microbiome Analyses ------------------------------------------------------


# Alpha-Diversity Analysis

library(nortest)

adiv <- gms %>% 
        estimate_richness(measures = c('Observed', 'Shannon', 'Simpson', 'chao')) %>% 
        rownames_to_column("id") %>% 
        mutate_at(vars(id), ~ str_replace(., fixed("."), "-"))

adiv_sum <- adiv %>% 
        left_join(., data.frame(sample_data(gms)) %>% rownames_to_column('id')) %>% 
        #filter(time == 'T1') %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), 
                     names_to = 'measures', values_to = 'values') %>%
        group_by(treatment, measures, time) %>% 
        nest() %>% 
        mutate(norm.test = purrr::map(data, function(x) (ad.test(x$values)))) %>% 
        mutate(pval = purrr::map_dbl(norm.test, function(x) x$p.value))
# Simpson for A and B treatment fail normality        

LM_adiv <- adiv %>% 
        left_join(., data.frame(sample_data(gms)) %>% rownames_to_column('id')) %>%
        mutate(Simpson = log(Simpson)) %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), 
                     names_to = 'measures', values_to = 'values') %>% 
        group_by(measures) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) 
                lmerTest::lmer(values ~ treatment*time + (1|study.id), data = x))) %>% 
        mutate(T1_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`))


# Beta Diversity Analysis 

# XN vs Ctr at Baseline
t1_ord <- ordinate(ps_baseline, method = 'PCoA', distance = 'bray')
# Plot ordinations
plot_ordination(ps_baseline, t1_ord, color='treatment', shape='treatment') + 
        geom_point(size = 3) + 
        ggtitle('XN vs Ctr at Baseline')
distmat <- phyloseq::distance(ps_baseline, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_baseline))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$treatment)
permutest(bdisp) # p = 0.03
#unequal beta dispersion btw intervention groups at baseline
adonis2(distmat~treatment, data = mdata) # p = 0.15

# XN vs Ctr at T2
t2_ord <- ordinate(ps_t2, method = 'PCoA', distance = 'jaccard')
# Plot ordinations
plot_ordination(ps_t2, t2_ord, color='treatment', shape='treatment') + 
        geom_point(size = 3) + 
        ggtitle('XN vs Ctr at Week 2')
distmat <- phyloseq::distance(ps_t2, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_t2))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$treatment)
permutest(bdisp) # p = 0.26
adonis2(distmat~treatment, data = mdata) # p = 0.23

# XN vs Ctr at T3
t3_ord <- ordinate(ps_t3, method = 'PCoA', distance = 'jaccard')
# Plot ordinations
plot_ordination(ps_t3, t3_ord, color='treatment', shape='treatment') + 
        geom_point(size = 3) + 
        ggtitle('XN vs Ctr at Week 4')
distmat <- phyloseq::distance(ps_t3, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_t3))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$treatment)
permutest(bdisp) # p = 0.07
adonis2(distmat~treatment, data = mdata) # p = 0.21

# XN vs Ctr at T4
t4_ord <- ordinate(ps_t4, method = 'PCoA', distance = 'jaccard')
# Plot ordinations
plot_ordination(ps_t4, t4_ord, color='treatment', shape='treatment') + 
        geom_point(size = 3) + 
        ggtitle('XN vs Ctr at Week 6')
distmat <- phyloseq::distance(ps_t4, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_t4))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$treatment)
permutest(bdisp) # p = 0.19
adonis2(distmat~treatment, data = mdata) # p = 0.54

# XN vs Ctr at T5
t5_ord <- ordinate(ps_t5, method = 'PCoA', distance = 'jaccard')
# Plot ordinations
plot_ordination(ps_t5, t5_ord, color='treatment', shape='treatment') + 
        geom_point(size = 3) + 
        ggtitle('XN vs Ctr at Week 8')
distmat <- phyloseq::distance(ps_t5, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_t5))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$treatment)
permutest(bdisp) # p = 0.35
adonis2(distmat~treatment, data = mdata) # p = 0.52

# Time by XN group
time_ord <- ordinate(ps_treat, method = 'PCoA', distance = 'jaccard')
# Plot ordinations
plot_ordination(ps_treat, time_ord, color = 'time') +
        geom_point(size = 3)+
        ggtitle('XN-treated group by timepoint')
distmat <- phyloseq::distance(ps_treat, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_treat))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$time)
permutest(bdisp) # p = 0.925
adonis2(distmat~time, data = mdata) # p = 1.0


# Negative Binomial Differential Abundance Analysis ----------------------------

library(glmmTMB)
library(emmeans)
library(multcomp)

micro_counts <- ps_counts %>% 
        #set.seed(711)
        rarefy_even_depth(rngseed = 711) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        rownames_to_column('id') %>% 
        dplyr::select(-c(ASV1024, ASV1520))

mdata <- sample_data(ps_counts) %>%
        data.frame() %>%
        mutate(study.id = str_replace(study.id, "XMaS", "")) %>% 
        rownames_to_column('id') 

fulldata <- left_join(mdata, micro_counts) 

# Test Normality of clr-transformed data
count_sum <- fulldata %>% 
        dplyr::select(-c(ASV1791, ASV647)) %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'value') %>% 
        group_by(ASV, treatment, time) %>% 
        nest() %>% 
        mutate(norm.test = purrr::map(data, function(x) (ad.test(x$value)))) %>% 
        mutate(pval = purrr::map_dbl(norm.test, function(x) x$p.value))

glmer.nb.data <- fulldata %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'value') %>% 
        group_by(ASV) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) glmmTMB(value ~ treatment*time + (1|study.id), family=nbinom2, data = x))) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,0,0,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,1,0,0,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(T3_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,1,0,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(T4_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,0,1,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(T5_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,0,0,1), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT2 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,1,0,0,0,1,0,0,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT3 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,1,0,0,0,1,0,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT4 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,0,1,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT5 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,0,1,0,0,0,1), nrow = 1)))[[1,6]][1])) %>%
        ungroup() %>% 
        mutate(T1_BvA_adj = p.adjust(T1_BvA, method = 'BH')) %>% 
        mutate(T2_BvA_adj = p.adjust(T2_BvA, method = 'BH')) %>%
        mutate(T3_BvA_adj = p.adjust(T3_BvA, method = 'BH')) %>%
        mutate(T4_BvA_adj = p.adjust(T4_BvA, method = 'BH')) %>%
        mutate(T5_BvA_adj = p.adjust(T5_BvA, method = 'BH')) %>%
        mutate(B_T1vT2_adj = p.adjust(B_T1vT2, method = 'BH')) %>%
        mutate(B_T1vT3_adj = p.adjust(B_T1vT3, method = 'BH')) %>%
        mutate(B_T1vT4_adj = p.adjust(B_T1vT4, method = 'BH')) %>%
        mutate(B_T1vT5_adj = p.adjust(B_T1vT5, method = 'BH'))


# Enterotype Cluster Analysis --------------------------------------------------

library(cluster)

ps_treat <- ps %>% 
  subset_samples(treatment == 'B')

# Jensen-Shannon ordination of XN-treated group 
treat.ord <- ordinate(ps_treat, method = 'PCoA', distance = 'jsd')

# Extract dataframe of jsd ordinations
ord_df <- plot_ordination(ps_treat, treat.ord, justDF = T)

# Calculate partitioning around medoids (PAM) clustering with 3 medoids
out <- pam(ord_df[, 1:2], 3, metric = "euclidean") 
silh_values <- pam(ord_df[, 1:2], 3, metric = "euclidean", keep.diss = T) 
# ave sil width = 0.51
plot(silh_values)

micro_counts <- ps_counts %>% 
        #set.seed(711)
        rarefy_even_depth(rngseed = 711) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        rownames_to_column('id')

# Extract clustering + add base_clus variable + add rarefied ASV counts
# For reference: clus 1 = Prevotella, clus 2 = Bacteroides, clus 3 = Ruminococcus
clusters <- out$clustering %>% 
        as.data.frame %>%
        rownames_to_column('id') %>% 
        dplyr::rename(clus = ".") %>% 
        dplyr::mutate_at(vars(clus), as.factor) %>% 
        left_join(ord_df %>% rownames_to_column('id'), .) %>% 
        mutate(base_clus = ifelse(study.id == "113", "1", clus)) %>% 
        #mutate(base_clus = case_when(study.id == "113" ~ "1",
                                     #.default = clus)) #%>%
        #mutate(across(base_clus, as_factor)) #%>% 
        dplyr::mutate_at(vars(base_clus), as.factor) %>% 
        dplyr::relocate(clus, base_clus) %>% 
        left_join(., micro_counts) %>%
        column_to_rownames('id')

#write_rds(clusters, "~/Projects/XMaSI/Data/clusters.rds")


# Beta Diversity of Enterotype*time

# Time + ET in XN group 
ettime_ord <- ordinate(ps_treat, method = 'PCoA', distance = 'jaccard')
# Plot ordinations
plot_ordination(ps_treat, ettime_ord, color = 'time') +
        geom_point(size = 3)+
        ggtitle('XN-treated group by timepoint')
distmat <- phyloseq::distance(ps_treat, method = 'jaccard')
#Set seed for reproducibility
set.seed(111)
#Extract meta data
mdata <-data.frame(sample_data(ps_treat)) %>% 
        rownames_to_column('id') %>% 
        left_join(clusters %>% rownames_to_column('id') %>% dplyr::select(id, base_clus))
#Verify variance is equal between groups, an assumption for PERMANOVA
bdisp <- betadisper(distmat, mdata$time)
permutest(bdisp) # p = 0.925
adonis2(distmat~base_clus*time, data = mdata) 
# p = 1.0


# Differential Abundance by Enterotype cluster

clus.nb.data <- clusters %>% 
        #simplify model with 3 time points instead of 5
        dplyr::filter(time %in% c('T1', 'T3', 'T5')) %>% 
        # Remove ASV's with too many zeros, not enough information for model
        dplyr::select(-c(ASV126, ASV647, ASV658, ASV1146, ASV1290, ASV1740, ASV1847, ASV1949, ASV2121)) %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'value') %>% 
        group_by(ASV) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) glmmTMB(value ~ base_clus*time + (1|study.id), family=nbinom2, data = x))) %>% 
        # test contrast for each cluster comparing taxa at close-out to taxa at baseline
        mutate(clus1_t5vt1 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,0,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(clus2_t5vt1 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,1,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(clus3_t5vt1 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,0,1), nrow = 1)))[[1,6]][1])) %>% 
        ungroup() %>% 
        mutate(clus1_adj = p.adjust(clus1_t5vt1, method = 'BH')) %>% 
        mutate(clus2_adj = p.adjust(clus2_t5vt1, method = 'BH')) %>%
        mutate(clus3_adj = p.adjust(clus3_t5vt1, method = 'BH'))


# Spearmans's correlation analysis between sample sites 

fullxns <- xns_f %>% 
        left_join(., xns_p) %>% 
        left_join(., xns_u) %>% 
        left_join(., metas) %>% 
        filter(treatment == "B")

# stool x plasma
cor.test(fullxns$f_xn, fullxns$p_xn, method='spearman')
#rho=0.50, p < 0.001
cor.test(fullxns$f_ixn, fullxns$p_ixn, method='spearman')
#rho=0.69, p < 0.001
cor.test(fullxns$f_8pn, fullxns$p_8pn, method='spearman')
#rho=0.52, p < 0.001
cor.test(fullxns$f_dxn, fullxns$p_dxn, method='spearman')
#rho=-0.02, p = 0.86
cor.test(fullxns$f_ddxn, fullxns$p_ddxn, method='spearman')
#rho= NA, p= NA

# stool x urine
cor.test(fullxns$f_xn, fullxns$u_xn, method='spearman')
#rho=0.57, p < 0.001
cor.test(fullxns$f_ixn, fullxns$u_ixn, method='spearman')
#rho=0.69, p < 0.001
cor.test(fullxns$f_8pn, fullxns$u_8pn, method='spearman')
#rho=0.79, p < 0.001
cor.test(fullxns$f_dxn, fullxns$u_dxn, method='spearman')
#rho=070., p < 0.001
cor.test(fullxns$f_ddxn, fullxns$u_ddxn, method='spearman')
#rho=0.40, p < 0.001


# Metabolite generation by Enterotype Cluster ----------------------------------

metab.by.clus <- clusters %>% 
      dplyr::filter(time != 'T1') %>% 
      mutate(across(.cols = starts_with('f_'), ~log(.x))) %>% 
      pivot_longer(starts_with("f_"), names_to="metab", values_to="value") %>% 
      group_by(metab) %>% 
      nest() %>% 
      mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ base_clus + (1|study.id), data = x))) %>% 
      mutate(clus1vsclus2 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0, -1, 0)))$`Pr(>F)`)) %>% 
      mutate(clus2vsclus3 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0, 1, -1)))$`Pr(>F)`)) %>% 
      mutate(clus1vsclus3 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0, 0, -1)))$`Pr(>F)`))


# Metabolite concentration by Enterotype Cluster -------------------------------

conc.by.et1 <- clusters %>% 
        # filter out baseline & by cluster 1 (Prevotella)
        dplyr::filter(time != 'T1' & base_clus == '1') %>% 
        rename(f_metab_sum = metab_sum) %>% 
        # change XN metabolite concentration from ng/g to ug/g
        mutate(across(.cols = starts_with('f_'), ~.x/1000)) %>%
        # log transform XN metabolite values
        mutate(across(.cols = starts_with('f_'), ~log(.x)))

LMM.conc.by.et1 <- conc.by.et1 %>% 
        pivot_longer(starts_with('f_'), names_to = 'metab', values_to = 'conc') %>% 
        group_by(metab) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ (1|study.id), data = x))) %>% 
        mutate(value = purrr::map_dbl(model, function(x) as.data.frame(summary(x)$coef[[1]])[[1]])) %>% 
        mutate(lower = purrr::map_dbl(model, function(x) as.data.frame(confint(x))[[3,1]])) %>% 
        mutate(upper = purrr::map_dbl(model, function(x) as.data.frame(confint(x))[[3,2]])) %>% 
        mutate(exp_val = purrr::map_dbl(value, function(x) exp(x))) %>% 
        mutate(exp_low = purrr::map_dbl(lower, function(x) exp(x))) %>% 
        mutate(exp_upp = purrr::map_dbl(upper, function(x) exp(x)))

conc.by.et2 <- clusters %>% 
        # filter out baseline & by cluster 2 (Bacteroides)
        dplyr::filter(time != 'T1' & base_clus == '2') %>% 
        rename(f_metab_sum = metab_sum) %>% 
        # change XN metabolite concentration from ng/g to ug/g
        mutate(across(.cols = starts_with('f_'), ~.x/1000)) %>%
        # log transform XN metabolite values
        mutate(across(.cols = starts_with('f_'), ~log(.x)))

LMM.conc.by.et2 <- conc.by.et2 %>% 
        pivot_longer(starts_with('f_'), names_to = 'metab', values_to = 'conc') %>% 
        group_by(metab) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ (1|study.id), data = x))) %>% 
        mutate(value = purrr::map_dbl(model, function(x) as.data.frame(summary(x)$coef[[1]])[[1]])) %>% 
        mutate(lower = purrr::map_dbl(model, function(x) as.data.frame(confint(x))[[3,1]])) %>% 
        mutate(upper = purrr::map_dbl(model, function(x) as.data.frame(confint(x))[[3,2]])) %>% 
        mutate(exp_val = purrr::map_dbl(value, function(x) exp(x))) %>% 
        mutate(exp_low = purrr::map_dbl(lower, function(x) exp(x))) %>% 
        mutate(exp_upp = purrr::map_dbl(upper, function(x) exp(x)))

conc.by.et3 <- clusters %>% 
        # filter out baseline & by cluster 3 (Ruminococcus)
        dplyr::filter(time != 'T1' & base_clus == '3') %>% 
        rename(f_metab_sum = metab_sum) %>% 
        # change XN metabolite concentration from ng/g to ug/g
        mutate(across(.cols = starts_with('f_'), ~.x/1000)) %>%
        # log transform XN metabolite values
        mutate(across(.cols = starts_with('f_'), ~log(.x)))

LMM.conc.by.et3 <- conc.by.et3 %>% 
        pivot_longer(starts_with('f_'), names_to = 'metab', values_to = 'conc') %>% 
        group_by(metab) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ (1|study.id), data = x))) %>% 
        mutate(value = purrr::map_dbl(model, function(x) as.data.frame(summary(x)$coef[[1]])[[1]])) %>% 
        mutate(lower = purrr::map_dbl(model, function(x) as.data.frame(confint(x))[[3,1]])) %>% 
        mutate(upper = purrr::map_dbl(model, function(x) as.data.frame(confint(x))[[3,2]])) %>% 
        mutate(exp_val = purrr::map_dbl(value, function(x) exp(x))) %>% 
        mutate(exp_low = purrr::map_dbl(lower, function(x) exp(x))) %>% 
        mutate(exp_upp = purrr::map_dbl(upper, function(x) exp(x)))


#Metabolite concentrations by participant characteristics ---------------------- 

fullxns <- xns_f %>% 
        left_join(., xns_p) %>% 
        left_join(., xns_u) %>% 
        left_join(., metas) %>% 
        filter(treatment == "B")

#8PN + sex
x8pn_sex_cors <- fullxns %>% 
        pivot_longer(cols = ends_with('8pn'), names_to = 'metab', values_to = 'conc') %>% 
        group_by(metab) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) lmerTest::lmer(conc ~ time + sex + (1|study.id), data = x)),
               t2_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[2,5]),
               t3_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[3,5]),
               t4_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[4,5]),
               t5_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[5,5]),
               sex_male = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[6,5])
        ) %>% 
        ungroup() %>% 
        mutate(t2_adj = p.adjust(t2_pval, method = 'BH'),
               t3_adj = p.adjust(t3_pval, method = 'BH'),
               t4_adj = p.adjust(t4_pval, method = 'BH'),
               t5_adj = p.adjust(t5_pval, method = 'BH'),
               sex_adj = p.adjust(sex_male, method = 'BH')
        )

#8PN + BMI
x8pn_BMI_cors <- fullxns %>% 
        pivot_longer(cols = ends_with('8pn'), names_to = 'metab', values_to = 'conc') %>% 
        group_by(metab) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) lmerTest::lmer(conc ~ time + BMI_cat + (1|study.id), data = x)),
               t2_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[2,5]),
               t3_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[3,5]),
               t4_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[4,5]),
               t5_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[5,5]),
               BMI_pval = map_dbl(mod, function(x) as.data.frame(summary(x)$coefficients)[6,5])
        ) %>% 
        ungroup() %>% 
        mutate(t2_adj = p.adjust(t2_pval, method = 'BH'),
               t3_adj = p.adjust(t3_pval, method = 'BH'),
               t4_adj = p.adjust(t4_pval, method = 'BH'),
               t5_adj = p.adjust(t5_pval, method = 'BH'),
               BMI_adj = p.adjust(BMI_pval, method = 'BH')
        )


# Multi-Omic sPLS Analysis (XN Metabolites) ------------------------------------

# Test normality for xn metabs

xns_sum <- xns_f %>% 
        left_join(., metas) %>% 
        pivot_longer(cols = c(f_xn, f_ixn, f_8pn, f_ddxn), values_to = 'conc', names_to = 'metabolite') %>% 
        group_by(treatment, time, metabolite) %>% 
        nest() %>% 
        mutate(norm.test = purrr::map(data, function(x) (ad.test(x$conc)))) %>% 
        mutate(pval = purrr::map_dbl(norm.test, function(x) x$p.value))


# Data Prep --------------------------------------------------------------------

ps_spls <- ps_counts %>% 
        subset_samples(treatment == "B") %>% 
        subset_samples(time != "T1")

micro_clr <- ps_spls %>% 
        #set.seed(711)
        rarefy_even_depth(rngseed = 711) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('id')

# Metadata
meta_full<- sample_data(ps_counts) %>%
        data.frame() %>%
        mutate(study.id = str_replace(study.id, "XMaS", "")) %>% 
        dplyr::select(study.id, treatment, time) %>% 
        rownames_to_column('id')

#Add XN metabolite data
full_block <- left_join(meta_full, micro_clr) %>% 
        left_join(xns_f) %>% 
        filter(treatment == "B") %>% 
        filter(time != "T1") %>% 
        dplyr::rename('XN' = 'f_xn',
                      'IXN' = 'f_ixn',
                      '6PN' = 'f_6pn',
                      '8PN' = 'f_8pn',
                      'DXN' = 'f_dxn',
                      'DDXN' = 'f_ddxn')

#Create the microbiome block
micro_block <- full_block %>% 
        dplyr::select(starts_with('ASV'), 'id') %>% 
        column_to_rownames('id')

#Create the metabolite block
metab_block <- full_block %>% 
        dplyr::select(c(XN, IXN, '6PN', '8PN',DXN, DDXN, id)) %>% 
        column_to_rownames('id')


# Model Design Tuning ----------------------------------------------------------

spls_xns <- spls(micro_block, metab_block, multilevel = as.factor(full_block$study.id),
                 ncomp = 2, mode = 'regression', scale = T)
plotVar(spls_xns, cutoff = 0.3, title = 'Microbiome vs XN Metabolites', legend = c('Microbiome', 'Metabolites'), 
        var.names = T, col = c('darkorchid', 'brown'))
cor(spls_xns$variates$X, spls_xns$variates$Y) #0.76, 0.85


#Tune the number of components
perf_xns <- perf(spls_xns, validation = 'Mfold', folds = 10, nrepeat = 5)
plot(perf_xns, criterion = 'Q2.total') # 1 Latent Component

# Tune variable parameters to keep
list.keepX <- c(seq(10, 50, 5))
list.keepY <- c(6) # include all XNs variables since there are only 6

tune.spls.xns <- tune.spls(X = micro_block, Y = metab_block, ncomp = 2,
                           multilevel = as.factor(full_block$study.id),
                           test.keepX = list.keepX,
                           test.keepY = list.keepY,
                           nrepeat = 5, folds = 10, 
                           mode = 'regression', measure = 'cor'
)
plot(tune.spls.xns)

optimal.keepX <- tune.spls.xns$choice.keepX # Micro: 50, 10
optimal.keepY <- tune.spls.xns$choice.keepY # Metabolites: 6, 6
optimal.ncomp <- length(optimal.keepX)


#Final Model -------------------------------------------------------------------

final_model <- spls(micro_block, metab_block,
                    multilevel = full_block$study.id, #accounts for repeated measurements
                    ncomp = optimal.ncomp,
                    keepX = optimal.keepX,
                    keepY = optimal.keepY,
                    mode = 'regression')
plotVar(final_model, cutoff = 0.3, title = 'Microbiome vs XN Metabolites', legend = c('Microbiome', 'Metabolites'), 
        var.names = T, col = c('darkorchid', 'brown'))
plotLoadings(final_model, comp = 1, ndisplay = 20, size.name = 1.3)
selectVar(final_model, comp = 1)
cor(final_model$variates$X, final_model$variates$Y) #0.78, 0.796 #with ASV-level 0.81, 0.79

final_ASV <- selectVar(final_model)$X$value %>% 
        #filter(abs(value.var) >= 0.1) %>% 
        rownames()


final_loadings_ASV <- data.frame(tax_table(ps)) %>% 
        rownames_to_column('ASV') %>% 
        right_join(selectVar(final_model)$X$value %>% rownames_to_column('ASV')) #%>% 
filter(abs(value.var) >= 0.1)



# Create heatmap
ASV_names <- final_model$X %>% 
        colnames() %>% 
        as.data.frame() %>% 
        rename(ASV = ".") %>% 
        right_join(., alltax) %>% 
        dplyr::select(Genus) %>% 
        as.vector() %>%
        unlist() %>% 
        as.character()

dev.off()
cim(final_model, comp = 1:2, 
    xlab = "Genera", ylab = 'Metabolites', 
    margins = c(14,9),
    scale = T,
    transpose=T,
    zoom = T, 
    row.names = ASV_names,
    cutoff = 0.25)


# Bile Acid Analysis -----------------------------------------------------------

LMM_ba <- metas %>% 
        left_join(bas) %>% 
        dplyr::mutate(across(starts_with('BA_'), ~log(.x + 1))) %>% 
        dplyr::select(treatment, time, study.id, BA_CP, BA_UCP, BA_AUC, BA_UC, BA_UCS) %>% 
        pivot_longer(starts_with("BA_"), names_to = 'bas', values_to = 'conc') %>% 
        group_by(bas) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ treatment*time + (1|study.id), data = x))) %>%
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>%
        ungroup() %>%
        mutate(t1_BvA_adj = p.adjust(T1_BvA, method = 'BH')) %>% 
        mutate(t2_BvA_adj = p.adjust(T2_BvA, method = 'BH')) %>%
        mutate(t3_BvA_adj = p.adjust(T3_BvA, method = 'BH')) %>%
        mutate(t4_BvA_adj = p.adjust(T4_BvA, method = 'BH')) %>%
        mutate(t5_BvA_adj = p.adjust(T5_BvA, method = 'BH')) %>%
        mutate(t2_BvB_adj = p.adjust(T2_BvB, method = 'BH')) %>%
        mutate(t3_BvB_adj = p.adjust(T3_BvB, method = 'BH')) %>%
        mutate(t4_BvB_adj = p.adjust(T4_BvB, method = 'BH')) %>%
        mutate(t5_BvB_adj = p.adjust(T5_BvB, method = 'BH'))


# Bile Acid by Enterotype ------------------------------------------------------

ba.by.base.clus <- clusters %>% 
        rownames_to_column('id') %>% 
        left_join(bas) %>% 
        dplyr::select(id, base_clus, time, study.id, BA_AUC, BA_UCS) %>% 
        dplyr::mutate(across(starts_with('BA_'), ~log(.x +  1))) %>% 
        pivot_longer(starts_with("BA_"), names_to = "metab", values_to="value") %>% 
        group_by(metab) %>% 
        nest() %>% 
        dplyr::mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ base_clus*time + (1|study.id), data = x))) %>% 
        mutate(E1_T2vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E1_T3vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E1_T4vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E1_T5vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E2_T2vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E2_T3vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E2_T4vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E2_T5vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,0,1,0,0,0,1,0,0,0,0)))$`Pr(>F)`)) %>%
        mutate(E3_T2vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(E3_T3vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>% 
        mutate(E3_T4vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(E3_T5vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        ungroup() %>%
        mutate(E1_T2_adj = p.adjust(E1_T2vT1, method = 'BH')) %>% 
        mutate(E1_T3_adj = p.adjust(E1_T3vT1, method = 'BH')) %>%
        mutate(E1_T4_adj = p.adjust(E1_T4vT1, method = 'BH')) %>%
        mutate(E1_T5_adj = p.adjust(E1_T5vT1, method = 'BH')) %>%
        mutate(E2_T2_adj = p.adjust(E2_T2vT1, method = 'BH')) %>%
        mutate(E2_T3_adj = p.adjust(E2_T3vT1, method = 'BH')) %>%
        mutate(E2_T4_adj = p.adjust(E2_T4vT1, method = 'BH')) %>%
        mutate(E2_T5_adj = p.adjust(E2_T5vT1, method = 'BH')) %>% 
        mutate(E3_T2_adj = p.adjust(E3_T2vT1, method = 'BH')) %>%
        mutate(E3_T3_adj = p.adjust(E3_T3vT1, method = 'BH')) %>%
        mutate(E3_T4_adj = p.adjust(E3_T4vT1, method = 'BH')) %>%
        mutate(E3_T5_adj = p.adjust(E3_T5vT1, method = 'BH'))

# New group with clusters 1+3 and 2 (exploratory analysis - looks like the response in BA conc is coming from Prev and Rumin clusters)
ba.by.clus <- clusters %>% 
        rownames_to_column('id') %>% 
        left_join(bas) %>% 
        dplyr::select(id, base_clus, time, study.id, BA_AUC, BA_UCS) %>% 
        # Create new cluster variable - 'a' = Prev & Rumin, 'b' = Bacter
        mutate(new_clus = ifelse(base_clus %in% c('1', '3'), 'a', 'b')) %>% 
        dplyr::mutate(across(starts_with('BA_'), ~log(.x +  1))) %>% 
        pivot_longer(starts_with("BA_"), names_to = "metab", values_to="value") %>% 
        group_by(metab) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ new_clus*time + (1|study.id), data = x))) %>% 
        mutate(Ea_T2vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(Ea_T3vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(Ea_T4vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(Ea_T5vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(Eb_T2vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(Eb_T3vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>% 
        mutate(Eb_T4vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(Eb_T5vT1 = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>%
        mutate(Ea_T2_adj = p.adjust(Ea_T2vT1, method = 'BH')) %>% 
        mutate(Ea_T3_adj = p.adjust(Ea_T3vT1, method = 'BH')) %>%
        mutate(Ea_T4_adj = p.adjust(Ea_T4vT1, method = 'BH')) %>%
        mutate(Ea_T5_adj = p.adjust(Ea_T5vT1, method = 'BH')) %>%
        mutate(Eb_T2_adj = p.adjust(Eb_T2vT1, method = 'BH')) %>%
        mutate(Eb_T3_adj = p.adjust(Eb_T3vT1, method = 'BH')) %>%
        mutate(Eb_T4_adj = p.adjust(Eb_T4vT1, method = 'BH')) %>%
        mutate(Eb_T5_adj = p.adjust(Eb_T5vT1, method = 'BH')) 


# Short-chain Fatty Acid Analysis ----------------------------------------------

LMM_scfa <- metas %>% 
        left_join(scfas) %>% 
        dplyr::mutate(across(starts_with('FA_'), ~log(.x + 0.001))) %>% 
        dplyr::filter(!study.id %in% c(120, 125, 128)) %>% 
        pivot_longer(starts_with("FA_"), names_to = 'scfas', values_to = 'conc') %>% 
        group_by(scfas) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ treatment*time + (1|study.id), data = x))) %>%
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>%
        ungroup() %>%
        mutate(t1_BvA_adj = p.adjust(T1_BvA, method = 'BH')) %>% 
        mutate(t2_BvA_adj = p.adjust(T2_BvA, method = 'BH')) %>%
        mutate(t3_BvA_adj = p.adjust(T3_BvA, method = 'BH')) %>%
        mutate(t4_BvA_adj = p.adjust(T4_BvA, method = 'BH')) %>%
        mutate(t5_BvA_adj = p.adjust(T5_BvA, method = 'BH')) %>%
        mutate(t2_BvB_adj = p.adjust(T2_BvB, method = 'BH')) %>%
        mutate(t3_BvB_adj = p.adjust(T3_BvB, method = 'BH')) %>%
        mutate(t4_BvB_adj = p.adjust(T4_BvB, method = 'BH')) %>%
        mutate(t5_BvB_adj = p.adjust(T5_BvB, method = 'BH'))


# Gut Health Markers Analysis --------------------------------------------------

LMM_gmarks <- metas %>% 
        left_join(gmarks) %>%
        dplyr::filter(!id %in% c('s117-5', 's126-3', 's126-4', 's126-5')) %>%
        dplyr::mutate(across(GM_IL12:GM_plasma_LBP, ~as.numeric(.))) %>% 
        dplyr::mutate(across(starts_with('GM_'), ~log(.x + 0.001))) %>% 
        pivot_longer(starts_with("GM_"), names_to = 'GM', values_to = 'conc') %>% 
        group_by(GM) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ treatment*time + (1|study.id), data = x))) %>%
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>%
        ungroup() %>%
        mutate(t1_BvA_adj = p.adjust(T1_BvA, method = 'BH')) %>% 
        mutate(t2_BvA_adj = p.adjust(T2_BvA, method = 'BH')) %>%
        mutate(t3_BvA_adj = p.adjust(T3_BvA, method = 'BH')) %>%
        mutate(t4_BvA_adj = p.adjust(T4_BvA, method = 'BH')) %>%
        mutate(t5_BvA_adj = p.adjust(T5_BvA, method = 'BH')) %>%
        mutate(t2_BvB_adj = p.adjust(T2_BvB, method = 'BH')) %>%
        mutate(t3_BvB_adj = p.adjust(T3_BvB, method = 'BH')) %>%
        mutate(t4_BvB_adj = p.adjust(T4_BvB, method = 'BH')) %>%
        mutate(t5_BvB_adj = p.adjust(T5_BvB, method = 'BH'))

