# Script to run main pgls models between flower maleness against 4 alternative 
# predictors of insect evolutionary pressure.
# Last run: 2020.10.06

# Load packages -----------------------------------------------------------
library(tidyverse) # CRAN v1.3.0 
library(sensiPhy)  # CRAN v0.8.5
library(phylolm)   # CRAN v2.6.2   
library(patchwork) # CRAN v1.0.1 
library(broom)     # CRAN v0.7.0     
library(writexl)   # CRAN v1.3   
library(ggtext)    # [github::wilkelab/ggtext] v0.1.0.9000    
library(ggExtra)   # CRAN v0.9   
library(psych)     # CRAN v2.0.7     
library(GGally)    # CRAN v2.0.0    

# source local functions--------------------------------------------------------
source("R/zzz_functions.R")

# Load data ---------------------------------------------------------------
set.seed(1234522)
d <- read_csv("data/processed/full_data.csv")
t <- read.tree(file = "data/processed/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# Fit PGLS----------------------------------------------------------------------
# 1. maleness ~ nspe------------------------------------------------------------
m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda", boot = 1000)
summary(m1)

# save model
writexl::write_xlsx(x = tidy(m1), path = "output/supp/supp_tab_pgls_maleness_vs_nspe.xls")
saveRDS(object = m1, file = 'output/temp/pgls_maleness_vs_nspe.RDs')

# 2. maleness ~ nfam------------------------------------------------------------
m2 <- phylolm(maleness ~ log2(nfam+1), data = d, phy = t, model = "lambda", boot = 1000)
summary(m2) 

# save model
writexl::write_xlsx(x = tidy(m2), path = "output/supp/supp_tab_pgls_maleness_vs_nfam.xls")
saveRDS(object = m2, file = 'output/temp/pgls_maleness_vs_nfam.RDs')

# 3. maleness ~ ngui----------------------------------------------------------------
m3 <- phylolm(maleness ~ log2(ngui+1), data = d, phy = t, model = "lambda", boot = 1000)
summary(m3) 

# save model
writexl::write_xlsx(x = tidy(m3), path = "output/supp/supp_tab_pgls_maleness_vs_ngui.xls")
saveRDS(object = m3, file = 'output/temp/pgls_maleness_vs_ngui.RDs')

# 4. maleness ~ shan----------------------------------------------------------------
m4 <- phylolm(maleness ~ shan, data = d, phy = t, model = "lambda", boot = 1000)
summary(m4) 

# save model
writexl::write_xlsx(x = tidy(m4), path = "output/supp/supp_tab_pgls_maleness_vs_shan.xls")
saveRDS(object = m4, file = 'output/temp/pgls_maleness_vs_shan.RDs')
