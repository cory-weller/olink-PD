#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggbeeswarm)

counts_files <- list.files(path='counts/COUNTS', pattern='*.csv', full.names=TRUE)
npx_files    <- list.files(path='counts/NPX',    pattern='*.csv', full.names=TRUE)

dat <- foreach(f=counts_files, .combine='rbind') %do% {
    d <- fread(f, header=TRUE)
    split_f <- strsplit(f, split='_')[[1]]
    sample_type <- split_f[4]
    panel <- split_f[5]
    d[, 'sample_type' := toupper(sample_type)]
    d[, 'panel' := toupper(panel)]
    return(d[])
}

dat <- dat[COUNT >= 500]

dat[, COUNT2 := COUNT/EXTENSIONCONTROLCOUNT]
dat[, LOG2COUNT := log2(COUNT2)]
dat[, BATCHID := .GRP, by = list(ASSAY,PLATEID,PANEL)]

dat <- merge(dat, 
             dat[PATNO %like% "PLATE_CTRL", list(.N, 'CTRLMED'=median(LOG2COUNT)), by=BATCHID],
             by.x='BATCHID',
             by.y='BATCHID')[N==3]



dat[PATNO %like% "NEG_CTRL", list(.N), by=BATCHID][N != 3]
dat[PATNO %like% "CONTROL_SAMPLE", list(.N), by=BATCHID][N != 2]
dat[, NPX := LOG2COUNT - CTRLMED]

dat[, grp := paste0(sample_type, '_', panel)]


dat <- dat[! PATNO %like% 'CTRL'][! PATNO %like% 'CONTROL']

global_med <- median(dat$NPX)

dat <- merge(dat,
             dat[, list('ASSAY_MED_NPX' = median(NPX)), by=list(ASSAY)],
             by.x='ASSAY',
             by.y='ASSAY')

dat[, NPX := NPX - ASSAY_MED_NPX + global_med ]

dat.counts <- dat

dat.counts[, 'TYPE' := 'COUNTS']


# get NPX tables

dat.npx <- foreach(f=npx_files, .combine='rbind') %do% {
    d <- fread(f, header=TRUE)
    split_f <- strsplit(f, split='_')[[1]]
    sample_type <- split_f[4]
    panel <- split_f[5]
    d[, 'sample_type' := toupper(sample_type)]
    d[, 'panel' := toupper(panel)]
    return(d[])
}

dat.npx[, 'TYPE' := 'NPX']

dat.all <- rbindlist(list(dat.counts, dat.npx))


ggplot(dat, aes(x=grp, y=intNormNPX)) + geom_violin()

dat.sub <- dat[, c('PATNO','PANEL','UNIPROT','sample_type','grp','intNormNPX')]
dat.sub[, id := .GRP, by=list(PATNO,PANEL,UNIPROT,sample_type)]

dat.sub <- unique(dat.sub[, list(PATNO, PANEL, UNIPROT, sample_type, NPX=mean(intNormNPX)), by=id])
dat.wide <- dcast(dat.sub, PATNO+PANEL+UNIPROT~sample_type, value.var='NPX')
setnames(dat.wide, 'PANEL', 'panel')
ggplot(dat.wide, aes(x=CSF, y=PLASMA)) + geom_point() + facet_grid(.~panel)

dat[, grpMedian := median(logCOUNT), by=grp]

plasma_median <- median(dat[sample_type=='PLASMA']$logCOUNT)
dat[sample_type=='PLASMA', shiftedCOUNT := plasma_median + logCOUNT - grpMedian]

csf_median <- median(dat[sample_type=='CSF']$logCOUNT)
dat[sample_type=='CSF', shiftedCOUNT := csf_median + logCOUNT - grpMedian]

ggplot(dat, aes(x=grp, y=shiftedCOUNT)) + geom_boxplot()

dat[, Z := (shiftedCOUNT-median(shiftedCOUNT)/mad(shiftedCOUNT)), by=grp]

ggplot(dat, aes(x=grp, y=Z)) + geom_boxplot()


robust_z_i = (x[i]-median(x))/mad(x)

dat[, 'extNPX' := log2(COUNT/EXTENSIONCONTROLCOUNT) , by=list(ASSAY,PLATEID)]

ggplot(dat, aes(x=grp, y=NPX_2)) + geom_boxplot()

dat[, NPX_1 := extNPX -median(extNPX), by=PLATEID]
dat[, NPX_2 := NPX_1 -median(NPX_1), by=grp]

dat[, medNPX_1 := median(NPX_1), by=PLATEID]
dat[, NPX_2 := NPX_1 -median(extNPX), by=PLATEID]

NPX_i, j = ExtNPX_i, j-median (ExtNPX(Plate Controls_i))

dat[PATNO %like% "PLATE_CTRL"]