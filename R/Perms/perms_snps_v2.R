# perms for SNP associations

library(qtl2)
set.seed(67179474)

# read stuff
do <- readRDS("../cache/do.rds")
pr <- readRDS("../cache/probs.rds")
k <- readRDS("../cache/kinship.rds")

# pseudomarker maps
gmap <- insert_pseudomarkers(do$gmap, stepwidth="max", step=0.2)
pmap <- interp_map(gmap, do$gmap, do$pmap)

# phenotypes and covariates
phe <- log10(do$pheno[,2])
covar <- cbind(sex=(do$is_female*1),
               wbc=log10(do$pheno[,1]))

# connection to snp database
qv <- create_variant_query_func("~/Data/CCdb/cc_variants.sqlite")

# autosomes
# snpinfo -> snp probs

# generate all of the permutations (stratify on sex)
n_perm <- 1000
perms <- qtl2:::gen_strat_perm(n_perm, rownames(do$geno[[1]]), do$is_female)

# mclapply...
library(parallel)
pr_copy <- pr
cores <- 12 # detectCores()
opermA <- mclapply(1:n_perm, function(i) {
    for(j in 1:19) rownames(pr_copy[[j]]) <- rownames(pr[[j]])[perms[,i]]
    out <- scan1snps(pr_copy, pmap, phe, k, addcovar=covar, query_func=qv, chr=1:19, cores=1)
    maxlod(out$lod) }, mc.cores=cores)

broman::done("perms_snps_v2: done with autosome")

# all that again for the X chromosome
# X chr
L <- chr_lengths(gmap, TRUE)
n_permX <- n_perm*L[1]/L[2]
perms <- qtl2:::gen_strat_perm(n_permX, rownames(do$geno[[1]]), do$is_female)
opermX <- mclapply(1:n_permX, function(i) {
    rownames(pr_copy[["X"]]) <- rownames(pr[["X"]])[perms[,i]]
    out <- scan1snps(pr_copy, pmap, phe, k, addcovar=covar, query_func=qv, chr="X", cores=1)
    maxlod(out$lod) }, mc.cores=cores)

result <- list(A=cbind(pheno1=unlist(opermA)), X=cbind(pheno1=unlist(opermX)))
class(result) <- c("scan1perm", "list")
attr(result, "chr_lengths") <- L

saveRDS(result, "operm_snps_v2.rds")

broman::done("done with perms_snps_v2")
