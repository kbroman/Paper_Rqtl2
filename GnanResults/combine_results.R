# combine the kover results into a single file

phe <- c("av.weight", "fruitlength", "ttl.seedsPfruit")
result <- NULL
for(p in phe) {
    file <- paste0(p, ".scan.txt")
    x <- data.table::fread(file, data.table=FALSE)
    x <- x[!is.na(x$from.bp),]

    stopifnot( all(x$from == x$to) )
    stopifnot( all(x$from.marker == x$to.marker) )
    stopifnot( all(x$from.bp == x$to.bp) )

    this <- x[,c("from.marker", "from.bp", "chromosome", "logP")]
    colnames(this) <- c("marker", "bp", "chr", "logP")

    if(is.null(result)) result <- this
    else {
        stopifnot( all(this[,1:3] == result[,1:3]) )
        result <- cbind(result, this$logP)
    }
    colnames(result)[ncol(result)] <- p
}

# save to RDS file
saveRDS(result, "gnan_results.rds")
