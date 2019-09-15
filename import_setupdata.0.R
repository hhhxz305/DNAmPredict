library(TCGA2STAT)

#modified merge function that returns transposed data
OMICSBIND2<-function (dat1, dat2) 
{
    if (is.null(dat1) || is.null(dat2)) {
        message("Empy data object")
        return(NULL)
    }
    if ((class(dat1) == "data.frame" || class(dat1) == "matrix") & 
        (class(dat2) == "data.frame" || class(dat2) == "matrix")) {
        t1 <- sapply(colnames(dat1), function(s) nchar(s))
        t2 <- sapply(colnames(dat2), function(s) nchar(s))
        dat1.tumor <- dat1
        dat2.tumor <- dat2
        if (sum(t1 > 12) == ncol(dat1)) {
            dat1.type <- sapply(colnames(dat1), function(s) unlist(strsplit(s, 
                "-"))[4])
            dat1.tumor <- dat1[, grep("^01", dat1.type)]
            colnames(dat1.tumor) <- substr(colnames(dat1.tumor), 
                1, 12)
        }
        if (sum(t2 > 12) == ncol(dat2)) {
            dat2.type <- sapply(colnames(dat2), function(s) unlist(strsplit(s, 
                "-"))[4])
            dat2.tumor <- dat2[, grep("^01", dat2.type)]
            colnames(dat2.tumor) <- substr(colnames(dat2.tumor), 
                1, 12)
        }
        matching <- intersect(colnames(dat1.tumor), colnames(dat2.tumor))
        if (length(matching) == 0) {
            message("No matched samples")
            return(NULL)
        }
        dat1.good <- dat1.tumor[, matching]
        dat2.good <- dat2.tumor[, matching]
        rownames(dat1.good) <- paste("d1.", rownames(dat1.good), 
            sep = "")
        rownames(dat2.good) <- paste("d2.", rownames(dat2.good), 
            sep = "")
        mdata <- t(rbind(dat1.good, dat2.good))
        return(list(merged.data = t(mdata)))
    }
}

#downloading SNP and clinical data
luad.snp <- getTCGA(disease="LUAD", data.type="CNA_SNP", type="450K",clinical=TRUE)
cesc.snp <- getTCGA(disease="CESC", data.type="CNA_SNP", type="450K",clinical=TRUE)
coad.snp <- getTCGA(disease="COAD", data.type="CNA_SNP", type="450K",clinical=TRUE)
kirc.snp <- getTCGA(disease="KIRC", data.type="CNA_SNP", type="450K",clinical=TRUE)
brca.snp <- getTCGA(disease="BRCA", data.type="CNA_SNP", type="450K",clinical=TRUE)

#downloading methylation data
luad.met <- getTCGA(disease="LUAD", data.type="Methylation", type="450K",clinical=T)
cesc.met <- getTCGA(disease="CESC", data.type="Methylation", type="450K",clinical=T)
coad.met <- getTCGA(disease="COAD", data.type="Methylation", type="450K",clinical=T)
kirc.met <- getTCGA(disease="KIRC", data.type="Methylation", type="450K",clinical=T)
brca.met <- getTCGA(disease="BRCA", data.type="Methylation", type="450K",clinical=T)

#saving methylation for later use in gini/gap
save(luad.met,file='luad_met_raw.rdata')
save(cesc.met,file='cesc_met_raw.rdata')
save(coad.met,file='coad_met_raw.rdata')
save(kirc.met,file='kirc_met_raw.rdata')
save(brca.met,file='brca_met_raw.rdata')

luad.snp1 <- OMICSBIND2(dat1 =luad.snp$dat, dat2 =t(luad.snp$clinical))
LUAD11 <- OMICSBIND2(dat1 =luad.met$dat, dat2 =luad.snp1$merged.data)

cesc.snp1 <- OMICSBIND2(dat1 =cesc.snp$dat, dat2 =t(cesc.snp$clinical))
CESC11 <- OMICSBIND2(dat1 =cesc.met$dat, dat2 =cesc.snp1$merged.data)


coad.snp1 <- OMICSBIND2(dat1 =coad.snp$dat, dat2 =t(coad.snp$clinical))
COAD11 <- OMICSBIND2(dat1 =coad.met$dat, dat2 =coad.snp1$merged.data)

kirc.snp1 <- OMICSBIND2(dat1 =kirc.snp$dat, dat2 =t(kirc.snp$clinical))
KIRC11 <- OMICSBIND2(dat1 =kirc.met$dat, dat2 =kirc.snp1$merged.data)

brca.snp1 <- OMICSBIND2(dat1 =brca.snp$dat, dat2 =t(brca.snp$clinical))
BRCA11 <- OMICSBIND2(dat1 =brca.met$dat, dat2 =brca.snp1$merged.data)

#writing combined data by cancer type for second part of analysis
write.csv(LUAD11$merged.data,'LUAD.csv')
write.csv(CESC11$merged.data,'CESC.csv')
write.csv(COAD11$merged.data,'COAD.csv')
write.csv(KIRC11$merged.data,'KIRC.csv')
write.csv(BRCA11$merged.data,'BRCA.csv')