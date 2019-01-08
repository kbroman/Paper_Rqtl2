# figure 1: reproduce Fig 5 from Gatti et al (2014)
#
# "We regressed log neutrophil counts on founder allele dosages at
# each marker using a kinship correction with sex and log white blood cell counts as covariates"

# load the data; 2nd phenotype is "NEUT"
library(qtl2)
set.seed(33003221)
file <- "cache/do.rds"
if(file.exists(file)) {
    do <- readRDS(file)
} else {
    do <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Gatti2014/do.zip")
    saveRDS(do, file)
}

# pseudomarker maps
file <- "cache/maps_n_phe.RData"
if(file.exists(file)) {
    load(file)
} else {
    gmap <- insert_pseudomarkers(do$gmap, stepwidth="max", step=0.2)
    pmap <- interp_map(gmap, do$gmap, do$pmap)

    # phenotypes and covariates
    phe <- log10(do$pheno[,2])
    covar <- cbind(sex=(do$is_female*1),
                   wbc=log10(do$pheno[,1]))

    save(gmap, pmap, phe, covar, file=file)
}

# calculate genotype probabilities
file <- "cache/probs.rds"
if(file.exists(file)) {
#    pr <- readRDS(file)
} else {
    pr <- calc_genoprob(do, gmap, error_prob=0.002, map_function="c-f", cores=0)
    saveRDS(pr, file)
}

# calculate allele dosages
file <- "cache/aprobs.rds"
if(file.exists(file)) {
#    apr <- readRDS(file)
} else {
    apr <- genoprob_to_alleleprob(pr, cores=0)
    saveRDS(apr, file)
}

# calculate kinship
file <- "cache/kinship.rds"
if(file.exists(file)) {
    k <- readRDS(file)
} else {
    k <- calc_kinship(apr, "loco", cores=0)
    saveRDS(k, file)
}

# genome scan with full model
file <- "cache/out_full.rds"
if(file.exists(file)) {
    out_full <- readRDS(file)
} else {
    out_full <- scan1(pr, phe, k, addcovar=covar, cores=0)
    saveRDS(out_full, file)
}

# genome scan with additive model
file <- "cache/out_add.rds"
if(file.exists(file)) {
    out_add <- readRDS(file)
} else {
    out_add <- scan1(apr, phe, k, addcovar=covar, cores=0)
    saveRDS(out_add, file)
}

# functions to query snps and genes
qv <- create_variant_query_func("~/Data/CCdb/cc_variants.sqlite")
qg <- create_gene_query_func("~/Data/CCdb/mouse_genes_mgi.sqlite")

# GWAS at SNPs

file <- "cache/out_snps.rds"
if(file.exists(file)) {
    out_snps <- readRDS(file)
} else {
    out_snps <- scan1snps(pr, pmap, phe, k,
                          addcovar=covar, query_func=qv, cores=0)
    saveRDS(out_snps, file)
}

# estimate coefficients on chr 1
file <- "cache/coef_c1.rds"
if(file.exists(file)) {
    co <- readRDS(file)
} else {
    co <- scan1coef(apr[,1], phe, k[1], addcovar=covar)
    co[,1:8] <- co[,1:8] - rowMeans(co[,1:8])
    saveRDS(co, file)
}

file <- "cache/blup_c1.rds"
if(file.exists(file)) {
    blup <- readRDS(file)
} else {
    blup <- scan1blup(apr[,1], phe, k[1], addcovar=covar, cores=0)
    saveRDS(blup, file)
}

# snp scan of a small region
li <- lod_int(out_add, pmap, chr=1)
file <- "cache/out_finemap_c1.rds"
if(file.exists(file)) {
    out_finemap <- readRDS(file)
} else {
    out_finemap <- scan1snps(pr, pmap, phe, k, addcovar=covar,
                             query_func=qv, chr=1, start=li[1], end=li[3],
                             keep_all_snps=TRUE, cores=0)
    saveRDS(out_finemap, file)
}
genes <- qg(chr=1, start=li[1], end=li[3])

# subset genes to 1/2
sub <- sort(unique(c(sample(nrow(genes), nrow(genes)/2),
                     match(c("Cxcr4", "Tmcc2"), genes$Name))))
genes <- genes[sub,]
gene_col <- rep("gray40", nrow(genes))
names(gene_col) <- genes$Name
gene_col[c("Cxcr4", "Tmcc2")] <- linecolor


# load permutation results
operm_full <- readRDS("Perms/operm_full.rds")
operm_add <- readRDS("Perms/operm_add.rds")
operm_snps <- readRDS("Perms/operm_snps.rds")

# calculate thresholds
thr_full <- summary(operm_full)
thr_add <- summary(operm_add)
thr_snps <- summary(operm_snps)

# ylim
ymx_full <- thr_full$A/thr_add$A*maxlod(out_add)*1.04
ymx_add <- maxlod(out_add)*1.04
ymx_snps <- thr_snps$A/thr_add$A*maxlod(out_add)*1.04

# make the plots
altcolor <- "green4"
linecolor <- "violetred"
panel_lab_adj <- c(0.12, 0.06)
panel_lab_cex <- 1.3

res <- 256
for(figtype in c("png", "eps")) {
if(figtype == "png") {
    png("../Figs/fig1.png", height=7.5*res, width=10*res, pointsize=14, res=res)
} else {
    postscript("../Figs/fig1.eps", paper="special", height=7.5, width=10, horizontal=FALSE, onefile=FALSE)
    panel_lab_adj[1] <- 0.10
}
layout(cbind(rep(1:3, each=4),
             c(4,4,4,5,5,6,6,6,7,7,7,7)))
par(mar=c(2.1, 4.1, 1.6, 1.1))

plot(out_full, pmap, xlab="", ylim=c(0, ymx_full), altcol=altcolor)
u <- par("usr")
endA <- xpos_scan1(pmap, thechr=19, thepos=max(pmap[[19]]))+25/2
segments(u[1], thr_full$A, endA, thr_full$A, col=linecolor, lty=2)
segments(endA, thr_full$X, u[2], thr_full$X, col=linecolor, lty=2)
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "A", font=2, xpd=TRUE, cex=panel_lab_cex)

plot(out_add, pmap, xlab="", ylim=c(0, ymx_add), altcol=altcolor)
segments(u[1], thr_add$A, endA, thr_add$A, col=linecolor, lty=2)
segments(endA, thr_add$X, u[2], thr_add$X, col=linecolor, lty=2)
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "B", font=2, xpd=TRUE, cex=panel_lab_cex)

plot(out_snps$lod, out_snps$snpinfo, altcol=altcolor, xlab="",
     ylim=c(0, ymx_snps))
segments(u[1], thr_snps$A, endA, thr_snps$A, col=linecolor, lty=2)
segments(endA, thr_snps$X, u[2], thr_snps$X, col=linecolor, lty=2)
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "C", font=2, xpd=TRUE, cex=panel_lab_cex)

par(mar=c(0.5,4.1,1.6,1.1))
ymx <- max(abs(blup[,1:8]))*1.04 * 1.6
mgp <- c(2.1, 0.3, 0)
if(figtype=="eps") mgp <- c(2.8, 0.3, 0)
plot_coefCC(blup, pmap, xaxt="n", ylim=c(-ymx, ymx), xlab="", mgp=mgp)
legend("topleft", ncol=4, lwd=2, col=CCcolors, legend=names(CCcolors), bg="gray92")
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "D", font=2, xpd=TRUE, cex=panel_lab_cex)
par(mar=c(3.1,4.1,0,1.1))
plot(out_add, pmap[1], xlab="", xaxt="n", mgp=mgp)
axis(side=1, at=pretty(par("usr")[1:2]), tick=FALSE, mgp=c(0, 0.2, 0))
title(xlab="Chr 1 position (Mbp)", mgp=c(1.8, 0, 0))

par(mar=c(0.5,4.1,1.6,1.1))
plot(out_finemap$lod, out_finemap$snpinfo, show=TRUE, drop=1.5,
     xaxt="n", xlab="")
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "E", font=2, xpd=TRUE, cex=panel_lab_cex)

par(mar=c(3.1,4.1,0,1.1))
plot_genes(genes, col=gene_col, xlab="", mgp=c(0,0.2,0),
           xlim=c(u[1], u[2]), xaxs="i")
title(xlab="Chr 1 position (Mbp)", mgp=c(1.8, 0, 0))
dev.off()
} # end loop over figure type
