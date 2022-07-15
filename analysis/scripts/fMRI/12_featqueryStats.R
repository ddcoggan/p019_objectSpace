# Title     : TODO
# Objective : TODO
# Created by: dave
# Created on: 3/22/22

# p008 exp2
# reads in timeseries.csv and plots response profiles to 6 different conditions across various brain regions
# produces text file with relevant stats

packages <- c('ggplot2', 'plyr', 'multcomp', 'broom', 'lsr', 'ez', 'plotrix', 'nlme', 'grid', 'Dict', 'stringr')
lapply(packages, require, character.only = TRUE)

sig_code <- function(x){
  y <- c()
  y[x > .05] = ''
  y[x <= .05] = '*'
  y[x <= .01] = '**'
  y[x <= .001] = '***'
  return(y)
}

voxelSizes = c(4,16,64) # data file only contains analyses for 64 voxels currently

# results directory
resdir <- 'analysis/results/fMRI/ROIresponses/topUp/noB0/doubleGamma/allSubjects'

data <- data.frame(read.csv('data/fMRI/group/featquery.csv', header = T, sep=',', strip.white=T)) # get data

# analysis 1: ANOVAs on each nml region without cat voxels
outFile = file.path(resdir, "NML_ANOVA.txt")
sink(outFile, append=F)
sink()
for (region in c('nml_lat_1','nml_lat_2','nml_med_1','nml_med_2')){
  theseData = droplevels(data[startsWith(as.character(data$roi), region) & endsWith(as.character(data$roi), 'noCatVoxels') & data$condition %in% c('body','face','house','object'),])
  theseData$region = NA
  theseData$region = substr(theseData$roi, 1, 9)
  dataAggregate = aggregate(PSC ~ subject + condition, data = theseData, FUN = mean)
  model <- ezANOVA(data = dataAggregate, dv = PSC, wid = subject, within = .(condition), detailed = TRUE)
  model$ANOVA$pes = model$ANOVA$SSn/(model$ANOVA$SSn + model$ANOVA$SSd) # add partial eta squared

  # create report file
  sink(outFile, append=T)
  cat(sprintf('### %s ###\n\n', region))
  cat(sprintf('### ANOVA ###\n', region))
  cat(capture.output(model), sep = '\n')
  sink()

  ## pw contrasts
  condsA = c('body','body','body','face','face','house')
  condsB = c('face','house','object','house','object','object')

  pwp_unc <- pairwise.t.test(dataAggregate$PSC, dataAggregate$condition, paired = T, p.adjust.method = "none")
  df <- data.frame(condA = condsA, condB = condsB, meanA=NA, meanB=NA, meandiff=NA, df = length(levels(dataAggregate$subject))-1, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
  for (i in 1:length(condsA)){
    A <- dataAggregate$PSC[dataAggregate$condition == condsA[i]]
    B <- dataAggregate$PSC[dataAggregate$condition == condsB[i]]
    df$meanA[i] <- mean(A)
    df$meanB[i] <- mean(B)
    df$meandiff[i] <- mean(A-B)
    test <- t.test(A,B,paired=T)
    df$t[i] <- test$statistic[[1]]
    df$p[i] <- pwp_unc$p.value[rownames(pwp_unc$p.value) == condsB[i], colnames(pwp_unc$p.value) == condsA[i]]
    df$d[i] <- cohensD(A,B,method='paired')
  }
  df$p.sig = sig_code(df$p)
  df$p.fdr <- p.adjust(df$p, method = 'holm')
  df$p.fdr.sig = sig_code(df$p.fdr)

  # create report in resdir
  sink(outFile, append = T)
  cat(sprintf('### post hocs ###\n'))
  cat(capture.output(df), sep = '\n')
  cat('\n\n')
  sink()
}


# analysis 2: ANOVAs on each category selective region region
outFile = file.path(resdir, "catANOVA.txt")
sink(outFile, append=F)
sink()
for (region in c('eba','fba','ofa','ffa','psts','opa','ppa','rsc')){
  theseData = droplevels(data[startsWith(as.character(data$roi), region) & data$condition %in% c('animate_spiky','animate_stubby','inanimate_spiky','inanimate_stubby'),])
  theseData$region = NA
  theseData$region = substr(theseData$roi, 1, nchar(region))
  theseData$animacy = NA
  theseData$aspectRatio = NA
  for (row in 1:nrow(theseData)){
    conds = strsplit(as.character(theseData$condition[row]), '_')
    theseData$animacy[row] = conds[[1]][1]
    theseData$aspectRatio[row] = strsplit(conds[[1]][2], '>')[[1]][1]
  }
  dataAggregate = aggregate(PSC ~ subject + animacy + aspectRatio + condition, data = theseData, FUN = mean)
  model <- ezANOVA(data = dataAggregate, dv = PSC, wid = subject, within = .(animacy, aspectRatio), detailed = TRUE)
  model$ANOVA$pes = model$ANOVA$SSn/(model$ANOVA$SSn + model$ANOVA$SSd) # add partial eta squared

  # create report file
  sink(outFile, append=T)
  cat(sprintf('### %s ###\n\n', region))
  cat(sprintf('### ANOVA ###\n', region))
  cat(capture.output(model), sep = '\n')
  sink()

  ## pw contrasts
  condsA = c('animate_spiky','animate_spiky','animate_spiky','animate_stubby','animate_stubby','inanimate_spiky')
  condsB = c('animate_stubby','inanimate_spiky','inanimate_stubby','inanimate_spiky','inanimate_stubby','inanimate_stubby')

  pwp_unc <- pairwise.t.test(dataAggregate$PSC, dataAggregate$condition, paired = T, p.adjust.method = "none")
  df <- data.frame(condA = condsA, condB = condsB, meanA=NA, meanB=NA, meandiff=NA, df = length(levels(dataAggregate$subject))-1, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
  for (i in 1:length(condsA)){
    A <- dataAggregate$PSC[dataAggregate$condition == condsA[i]]
    B <- dataAggregate$PSC[dataAggregate$condition == condsB[i]]
    df$meanA[i] <- mean(A)
    df$meanB[i] <- mean(B)
    df$meandiff[i] <- mean(A-B)
    test <- t.test(A,B,paired=T)
    df$t[i] <- test$statistic[[1]]
    df$p[i] <- pwp_unc$p.value[rownames(pwp_unc$p.value) == condsB[i], colnames(pwp_unc$p.value) == condsA[i]]
    df$d[i] <- cohensD(A,B,method='paired')
  }
  df$p.sig = sig_code(df$p)
  df$p.fdr <- p.adjust(df$p, method = 'holm')
  df$p.fdr.sig = sig_code(df$p.fdr)

  # create report in resdir
  sink(outFile, append = T)
  cat(sprintf('### post hocs ###\n'))
  cat(capture.output(df), sep = '\n')
  cat('\n\n')
  sink()
}
