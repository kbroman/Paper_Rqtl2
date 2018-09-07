# perms for additive allele model

library(qtl2)
set.seed(12435868)

# read stuff
do <- readRDS("../cache/do.rds")
apr <- readRDS("../cache/aprobs.rds")
k <- readRDS("../cache/kinship.rds")

# pseudomarker maps
gmap <- insert_pseudomarkers(do$gmap, stepwidth="max", step=0.2)
pmap <- interp_map(gmap, do$gmap, do$pmap)

# phenotypes and covariates
phe <- log10(do$pheno[,2])
covar <- cbind(sex=(do$is_female*1),
               wbc=log10(do$pheno[,1]))

n_perm <- 1000
operm_add <- scan1perm(apr, phe, k, addcovar=covar, cores=0, n_perm=n_perm, perm_Xsp=TRUE,
                       chr_lengths = chr_lengths(gmap))
saveRDS(operm_add, "operm_add.rds")
