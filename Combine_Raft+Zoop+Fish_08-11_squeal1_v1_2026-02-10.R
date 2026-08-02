# Combine Covariates for Squeal 1 logistic regressions
# 29 January 2026

rm(list = ls()) 
graphics.off()

# load data from 
# /Squeal1_Spline_dailies_DDJ+EPF/Read+Spline_Chl_MARSS-DLM_result_27Jan2026.R
# Langevin was fitted to smoothed predix daily chl from MARSS-DLM (x0 in mat1)
#save(dat0c,mat1,Stut,dat4reg,Tindx,file=Fname)
load(file='Spline_MARSSdlm_predix-tprobs_squeal1.Rdata')
#
print('covariates for 2008-2011 gathered to 27 January 2026',quote=F)
print(dim(dat0c))
print(dat0c[1,])

#print('From spline: distance to threshold stats',quote=F)
#print(dim(Stut))
#print(Stut[1,])

# load Daphnia size data from Cascade_Archive\Squeal-zoops-daily
# Dbar is Daphnia mean length mm, Dvar is variance, DCV is CV
#save(ZDat,file='Peter_Daphnia_length_daily2008-2011.Rdata')
load(file='Peter_Daphnia_length_daily2008-2011.Rdata')
print('Daphnia length data',quote=F)
print(dim(ZDat))
print(ZDat[1,])

# Load LMB + pelagic minnow index from 2008-2011
# source Build_Fish_Dataset_2008-2011_2026-02-10.R 
#   built from Build_Minnow_Index_2008-2011_28Jan2026.R
#   and Build_LMB_data_Peter_2008-2011.R
# save(fish0,file='LMB+minnow_2008-2011.Rdata')
load(file='LMB+minnow_2008-2011.Rdata')
print('minnow and lmb ratio')
print(fish0[1,])
print(dim(fish0))

# begin merge of covariates ===================================================

# merge dat0c with Daphnia length
dat1c = merge(x=dat0c,y=ZDat,by.x=c('year','doy'),all.x=T)
print('merge of initial covariates with Daphnia stats',quote=F)
print(dim(dat1c))
print(dim(na.omit(dat1c)))
print(dat1c[1,])

# merge dat1c (all limno + zoop covariates) with fish data
dat2c = merge(x=dat1c,y=fish0,by=c('year','doy'))
print('merge of dat2c with minnow + lmb data',quote=F)
print(dim(dat2c))
print(dim(na.omit(dat2c)))
print(dat2c[1, ])
print('totminnow are NA in 2011 due to missing pelagic data',quote=F)
print('the missing data are probably zeroes',quote=F)
print('NA totminnow were replace with 0 as placeholder',quote=F)
print('they should be replaced with litmin after merge',quote=F)
dat2c$totminnow = ifelse(dat2c$totminnow != 0,dat2c$totminnow,dat2c$litmin)
print('',quote=F)
print('repeating NA test for dat2c',quote=F)
print(dim(dat2c))
print(dim(na.omit(dat2c)))
print(dat2c[1, ])

# allcov is a dataframe with all covariates and no missing data
allcov = na.omit(dat2c)

# save daily covariate data from squeal 1
# allcov is all covariates with missing data removed
# dat0c is raft + ZB, dat1c includes Daphnia from Zdat, 
#    dat2c includes minnow + LMB data from fish0
print('saving all covariates for use with any squeal 1 transition series',quote=F)
save(allcov,dat0c,dat1c,dat2c,ZDat,fish0,file='covariates_2008-2011.Rdata')

# merge spline results with covariates
Stut2 = as.data.frame(cbind(mat1,Stut,Tindx))
print('attach year doy and Tindex to spline results',quote=F)
print(dim(Stut2))
print(Stut2[1,])

# merge Daphnia covariates with spline result
dat3 = merge(x=allcov,y=Stut2,by.x=c('year','doy'))
print('all covariates incl. Limno + Daphnia + fish + spine_result',quote=F)
print(dim(dat3))
print(dim(na.omit(dat3)))
print(dat3[1,])

all4reg = na.omit(dat3)

print('saving all covariates + t-stats for Squeal 1 logistic regression',quote=F)
save(all4reg,file='CombineCovarSplineSqueal1.Rdata')

