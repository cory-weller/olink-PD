#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)


dd <- '/data/CARD/projects/proteomicsXprogression_PD/data-allreleases/amp-pd-v3-proteomics/'
filenames <- list.files(path=dd, pattern="D0[12]_matrix_.*csv$", recursive=T, full.names=T)

# combine data sets and transform to long format
dat <- foreach(filename=filenames, .combine='rbind') %do% {
        mode <- strsplit(basename(filename), split='_matrix_|\\.csv')[[1]][2]
        dat.tmp <- fread(filename)
        dat.tmp.long <- melt(dat.tmp, 
                             measure.vars=colnames(dat.tmp)[2:ncol(dat.tmp)],
                             variable.name='sample')
        dat.tmp.long[, 'mode' := mode]
        return(dat.tmp.long[])
        }

# split 'variable' column into its separate values
dat[, c('dataset', 'participant_id', 'visit', 'sampletype', 'ppea', 'D') :=
        tstrsplit(sample, split='-')]

# remove redundant 'sample' column
dat[, sample := NULL]

# reduce file size

## is the number after BLM and SVM always month?
## What is difference between BLM and SVM?
## BL = Baseline?

dat[visit=="BLM0T1",  month:= 0 ]
dat[visit=="SVM3T1",  month:= 3 ]
dat[visit=="SVM6T1",  month:= 6 ]
dat[visit=="SVM9T1",  month:= 9 ]
dat[visit=="SVM12T1", month:= 12]
dat[visit=="SVM18T1", month:= 18]
dat[visit=="SVM24T1", month:= 24]
dat[visit=="SVM30T1", month:= 30]
dat[visit=="SVM36T1", month:= 36]
dat[visit=="SVM42T1", month:= 42]
dat[visit=="SVM48T1", month:= 48]
dat[visit=="SVM54T1", month:= 54]
dat[visit=="SVM60T1", month:= 60]
dat[visit=="SVM72T1", month:= 72]
dat[visit=="SVM84T1", month:= 84]
dat[visit=="SVM96T1", month:= 96]
dat[, visit := NULL]

# Simplify D02 to 2, D01 to 1
dat[D=='D02', D := '2']
dat[D=='D01', D := '1']

# Renumber participant_ids for file size reduction? 413 participants
# id_key <- data.table(participant_id=unique(dat$participant_id))
# id_key[, idx := 1:.N]
# setkey(id_key, idx)
# dat <- merge(dat, id_key, by.x='participant_id', by.y='participant_id')

fwrite(dat, file='reports/protein_quantification.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# # calculate composite for those with info at all three months
# hq_participants <- dat[! is.na(summary_score), .N, by=list(participant_id, visit_month)][N==3, participant_id]
# dat.hq <- dat[participant_id %in% hq_participants]


# dat.hq.total <- dat.hq[, list('total'=sum(summary_score)), by=list(participant_id, visit_month)]
# dat.hq.total[, jittered := total + runif(.N, -0.5, 0.5)]

# lowstart_pids <- unique(dat.hq.total[visit_month <= 6 & total <= 25, participant_id])

# ggplot(dat.hq.total[participant_id %in% lowstart_pids],
#     aes(x=visit_month, y=jittered, group=participant_id)) +
#     geom_line(alpha=0.1) +
#     labs(x='month', y='MDS composite score', title='score over time for initial score of 0')


###

# The following files relate to olink protein expression measured for participants


