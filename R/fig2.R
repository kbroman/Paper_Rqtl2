# Analysis of 19-way Arabidopsis MAGIC lines from Gnan et al. (2016)

library(qtl2)
library(broman)

# load the data
rdsfile <- "cache/arabmagic.rds"
if(file.exists(rdsfile)) {
    arabmagic <- readRDS(rdsfile)
} else {
    zipfile <- paste0("https://raw.githubusercontent.com/rqtl/",
                      "qtl2data/master/ArabMAGIC/arabmagic.zip")
    arabmagic <- read_cross2(zipfile)
    saveRDS(arabmagic, rdsfile)
}

# load the Gnan et al. results
gnan <- readRDS("../GnanResults/gnan_results.rds")

# subset the Gnan et al. results and arabmagic data to common markers
arabmagic <- pull_markers(arabmagic, gnan$marker %win% marker_names(arabmagic))
gnan <- gnan[gnan$marker %in% marker_names(arabmagic),]
# convert -logP to LOD scores
gnan_lod <- qchisq(-as.matrix(gnan[,4:6])*log(10), 18, lower.tail=FALSE, log.p=TRUE)/2/log(10)
rownames(gnan_lod) <- gnan$marker
class(gnan_lod) <- c("scan1", "matrix")

file <- "cache/arabmagic_lod.rds"
if(file.exists(file)) {
    out <- readRDS(file)
} else {
    pr <- calc_genoprob(arabmagic, error_prob=0.002, cores=0)
    out <- scan1(pr, arabmagic$pheno[,c("seed_weight", "fruit_length", "ttl_seedspfruit")], cores=0)
    saveRDS(out, file)
}

lab <- c("seed_size", "fruit_length", "seed_number_per_fruit")
tab5 <- lapply(lab, function(phe) data.table::fread(paste0("../GnanResults/table5_", phe, ".txt"), data.table=FALSE))
names(tab5) <- lab

# seed weight        1 @ 21.669 (row 2)
# fruit length       2 @ 11.207 (row 1)
# seed number/fruit  4 @  7.177 (row 7)
# - effects from Gnan et al
# - effects from R/qtl2
# - blup estimates from R/qtl2

# Estimates from Gnan et al. (2016), table 5
gnan_ests <- rbind(tab5[[1]][2,],
                   tab5[[2]][1,],
                   tab5[[3]][7,])[,-(1:2)]

# peak markers
mar <- c(rownames(max(gnan_lod, arabmagic$pmap, 1, 1)),
         rownames(max(gnan_lod, arabmagic$pmap, 2, 2)),
         rownames(max(gnan_lod, arabmagic$pmap, 3, 4)))

file <- "cache/arabmagic_ests.RData"
if(file.exists(file)) {
    load(file)
} else {
    pr <- calc_genoprob(arabmagic, error_prob=0.002, cores=0)

    # fit1 (plain estimates)
    out_fit1 <- lapply(1:3, function(i) fit1(pull_genoprobpos(pr, mar[i]), arabmagic$pheno[,c(2,8,4)[i],drop=FALSE]))
    fit1_ests <- matrix(unlist(sapply(out_fit1, "[", "coef")), byrow=TRUE, nrow=3)
    fit1_se <- matrix(unlist(sapply(out_fit1, "[", "SE")), byrow=TRUE, nrow=3)
    dimnames(fit1_ests) <- dimnames(fit1_se) <- list(colnames(arabmagic$pheno)[c(2,8,4)],
                                                 colnames(gnan_ests))
    fit1_lo <- fit1_ests - fit1_se*2
    fit1_hi <- fit1_ests + fit1_se*2

    # blup estimates
    out_blup <- lapply(1:3, function(i) scan1blup(list("1"=pr[[ c(1,2,4)[i] ]][,,mar[i],drop=FALSE]),
                                                  arabmagic$pheno[,c(2,8,4)[i], drop=FALSE], se=TRUE))

    blup_ests <- t(sapply(out_blup, function(a) a[1:19]+a[20]))
    blup_se <- t(sapply(out_blup, function(a) { b <- attr(a, "SE"); sqrt(b[1:19]^2+b[20]^2)}))
    dimnames(blup_ests) <- dimnames(blup_se) <- dimnames(fit1_ests)
    blup_lo <- blup_ests - blup_se*2
    blup_hi <- blup_ests + blup_se*2

    save(fit1_ests, fit1_se, fit1_lo, fit1_hi,
         blup_ests, blup_se, blup_lo, blup_hi, file=file)
}

yl_eff <- list(c(18.1, 25.8),
               c(11.7, 16.2),
               c(41.6, 59.9))

odd <- seq(1, 19, by=2)
even <- seq(2,19, by=2)

res <- 256
for(figtype in c("pdf", "eps")) {
if(figtype=="pdf")  {
    pdf("../Figs/fig2.pdf", height=7.5, width=10, pointsize=14)
} else {
    postscript("../Figs/fig2.eps", height=7.5, width=10, paper="special", horizontal=FALSE, onefile=FALSE)
}
par(mar=c(3.8, 3.5, 1.6, 1.1))
panel_lab_adj <- c(0.10, 0.06)
panel_lab_cex <- 1.3

if(figtype=="eps") panel_lab_adj[1] <- 0.08

blue <- "slateblue"
red <- "violetred"
green <- broman::brocolors("web")["green"]


par(mfcol=c(3,2))
yl <- c(0, max(maxlod(out), maxlod(gnan_lod))*1.05)
labels <- c("seed weight", "fruit length", "no. seeds per fruit")
for(i in 1:3) {
    plot(gnan_lod, arabmagic$gmap, gap=2, ylim=yl, main=labels[i], lod=i, xlab="", lwd=1, col=blue, mgp.y=c(1.8,0.3,0))
    plot(out, arabmagic$gmap, gap=2, add=TRUE, lod=i, lwd=1, col=red)
    u <- par("usr")
    text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], LETTERS[i], font=2, xpd=TRUE, cex=panel_lab_cex)

    if(i==1) {
        legend("topright", lwd=2, col=c(blue, red), c("Gnan et al.", "R/qtl2"), bg="gray92")
    }
}




# add chr and position to labels
marpos <- list(max(gnan_lod, arabmagic$pmap, 1, 1),
               max(gnan_lod, arabmagic$pmap, 2, 2),
               max(gnan_lod, arabmagic$pmap, 3, 4))
chr <- sapply(marpos, "[[", 1)
mbp <- broman::myround(sapply(marpos, "[[", 2), 1)

labels_wpos <- paste0(labels, " (chr ", chr, " at ", mbp, " Mbp)")



xd <- 0.2
for(i in 1:3) {
    grayplot(1:19-xd, gnan_ests[i,], ylim=yl_eff[[i]], xlab="", xat=NA, main=labels_wpos[i],
             ylab=paste("Mean", labels[i]), xlim=c(0.5, 19.5), xaxs="i",
             vlines=1:19, vlines.col="gray85", vlines.lwd=8, bg=blue, mgp.y=c(1.8,0.3,0))

    axis(side=1, at=odd, colnames(gnan_ests)[odd], tick=FALSE, mgp=c(0, 0.3, 0), las=2)
    axis(side=1, at=even, colnames(gnan_ests)[even], tick=FALSE, mgp=c(0, 0.3, 0), las=2)

    segments(1:19, fit1_lo[i,], 1:19,  fit1_hi[i,], lwd=1, col=red)
    points(1:19, fit1_ests[i,], pch=21, bg=red)

    segments(1:19+xd, blup_lo[i,], 1:19+xd,  blup_hi[i,], lwd=1, col=green)
    points(1:19+xd, blup_ests[i,], pch=21, bg=green)

    u <- par("usr")
    text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], LETTERS[i+3], font=2, xpd=TRUE, cex=panel_lab_cex)

    if(i==1) {
        legend("bottomleft", pch=21, pt.cex=1.3, pt.bg=c(blue, red, green),
               c("Gnan et al.", "R/qtl2 coef", "R/qtl2 blup"), bg="gray92")
    }
}

dev.off()
} # end loop over fig type
