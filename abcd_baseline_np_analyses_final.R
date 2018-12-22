#############################################################################################################
#############################################################################################################
## Code for:                                                                                               ##
## Thompson, W. K., Barch, D., Bjork, J., Gonzalez, R., Nagel, B., Nixon, S. J., & Luciana, M. (2018).     ## 
## The Structure of Cognition in 9 and 10 year-old Children and Associations with Problem Behaviors:       ##
## Findings from the ABCD Study’s Baseline Neurocognitive Battery. Developmental Cognitive Neuroscience.   ##
##                                                                                                         ##
## Wes Thompson (wes.stat@gmail.com)                                                                       ##
## Dec 22, 2018                                                                                            ##
#############################################################################################################
#############################################################################################################

####################
####################
## Load libraries ##
####################
####################

library(mvtnorm)
library(tableone)
library(parallel)
library(rstan)
library(loo)
library(gamm4)
library(Hmisc)
library(FactoMineR)
library(nFactors)
library(reshape2)
library(psych)
library(data.table)
library(mice)
library(abind)
library(cvTools)
library(modEvA)

####################################
####################################
## Load and manipulate nda17 data ##
####################################
####################################


## Read Rds file from DEAP (ABCD NDA version 1.1 release)
## scripts for creating the file "nda17_release1.1.Rds" can be found at: https://github.com/ABCD-STUDY/analysis-nda17

nda17 = readRDS("data/nda17_release1.1.Rds")

#########################################
#########################################
## Select & process data for analyses  ##
#########################################
#########################################

## Select demographics
ind_demog = c(which(names(nda17)=="age"),which(names(nda17)=="female"),which(names(nda17)=="race.ethnicity"),
			which(names(nda17)=="high.educ"),which(names(nda17)=="married"),which(names(nda17)=="household.income"))
names(nda17)[ind_demog]
summary(nda17[,ind_demog])

## Select nesting variables
ind_nest = c(which(names(nda17)=="abcd_site"),which(names(nda17)=="rel_family_id"));summary(nda17[,ind_nest])

nda17$abcd_site = as.character(nda17$abcd_site)
nda17$abcd_site[nda17$abcd_site=="site22"] = "site07"
nda17$abcd_site = factor(nda17$abcd_site)

## Select neuropsychological measures
ind_pea_ravlt = c(which(names(nda17)=="pea_ravlt_sd_trial_i_tc"),which(names(nda17)=="pea_ravlt_sd_trial_ii_tc"),
	which(names(nda17)=="pea_ravlt_sd_trial_iii_tc"),which(names(nda17)=="pea_ravlt_sd_trial_iv_tc"),
	which(names(nda17)=="pea_ravlt_sd_trial_v_tc")); names(nda17)[ind_pea_ravlt]; summary(nda17[,ind_pea_ravlt])
nda17$pea_ravlt_ld = apply(nda17[,ind_pea_ravlt],1,sum)

par(mfrow=c(1,2))
hist(nda17$pea_ravlt_ld)
hist(nda17$lmt_scr_perc_correct)

ind_np = c(which(names(nda17)=="nihtbx_picvocab_uncorrected"),which(names(nda17)=="nihtbx_flanker_uncorrected"),
	which(names(nda17)=="nihtbx_list_uncorrected"),which(names(nda17)=="nihtbx_cardsort_uncorrected"),
	which(names(nda17)=="nihtbx_pattern_uncorrected"),which(names(nda17)=="nihtbx_picture_uncorrected"),
	which(names(nda17)=="nihtbx_reading_uncorrected"),which(names(nda17)=="pea_ravlt_ld"),
	which(names(nda17)=="lmt_scr_perc_correct")); names(nda17)[ind_np]

## Select dependent measures 
ind_dv = c(which(names(nda17)=="cbcl_scr_syn_external_r"),
			which(names(nda17)=="cbcl_scr_syn_internal_r"),which(names(nda17)=="cbcl_scr_07_stress_r"))
names(nda17)[ind_dv]

######################
######################
## Rename variables ##
######################
######################

names(nda17)[which(names(nda17)=="age")] = "Age"; nda17$Age = nda17$Age/12
names(nda17)[which(names(nda17)=="female")] = "Female"
names(nda17)[which(names(nda17)=="race.ethnicity")] = "RaceEthnicity"
names(nda17)[which(names(nda17)=="high.educ")] = "HighestParentalEducation"
names(nda17)[which(names(nda17)=="married")] = "HouseholdMaritalStatus"
names(nda17)[which(names(nda17)=="household.income")] = "HouseholdIncome"
names(nda17)[which(names(nda17)=="abcd_site")] = "Site"
names(nda17)[which(names(nda17)=="rel_relationship")] = "Relationship"

names(nda17)[which(names(nda17)=="nihtbx_picvocab_uncorrected")] = "PicVocab"
names(nda17)[which(names(nda17)=="nihtbx_flanker_uncorrected")] = "Flanker"
names(nda17)[which(names(nda17)=="nihtbx_list_uncorrected")] = "List"
names(nda17)[which(names(nda17)=="nihtbx_cardsort_uncorrected")] = "CardSort"
names(nda17)[which(names(nda17)=="nihtbx_pattern_uncorrected")] = "Pattern"
names(nda17)[which(names(nda17)=="nihtbx_picture_uncorrected")] = "Picture"
names(nda17)[which(names(nda17)=="nihtbx_reading_uncorrected")] = "Reading"
names(nda17)[which(names(nda17)=="pea_ravlt_ld")] = "RAVLT"
names(nda17)[which(names(nda17)=="lmt_scr_perc_correct")] = "LMT"
names(nda17)[which(names(nda17)=="pea_wiscv_tss")] = "WISC-V"
names(nda17)[which(names(nda17)=="cbcl_scr_syn_external_r")] = "Externalizing"
names(nda17)[which(names(nda17)=="cbcl_scr_syn_internal_r")] = "Internalizing"
names(nda17)[which(names(nda17)=="cbcl_scr_07_stress_r")] = "Stress"

nda17$compl = "Incomplete"
nda17$compl[complete.cases(nda17[,c(ind_nest,ind_np)])] = "Complete"
table(nda17$compl)

## Create Table 1
vars1 <- c("Age", "Female", "RaceEthnicity", "HighestParentalEducation", "HouseholdMaritalStatus",
            "HouseholdIncome","Relationship","Site")
tab1 <- CreateTableOne(vars = vars1, data = nda17, strata = "compl")
tabAsStringMatrix <- print(tab1, printToggle = FALSE, noSpaces = TRUE)
tab1 = knitr::kable(tabAsStringMatrix)

## Create Table 2
vars2 <- c("PicVocab", "Flanker", "List", "CardSort", "Pattern",
            "Picture","Reading","RAVLT","WISC-V",
            "LMT","Externalizing","Internalizing","Stress")
tab2 <- CreateTableOne(vars = vars2, data = nda17,strata = "compl")
tabAsStringMatrix <- print(tab2, printToggle = FALSE, noSpaces = TRUE)
tab2 = knitr::kable(tabAsStringMatrix)

## Create SM Figures 1 & 2
pdf("figures_and_tables/fig_sm1.pdf")
par(mfrow=c(3,3))
for(p in 1:9){
	hist(scale(nda17[,ind_np[p]]), xlab="", freq=FALSE,main=names(nda17)[ind_np][p])
}
dev.off()

pdf("figures_and_tables/fig_sm2.pdf")
par(mfrow=c(2,2))
for(p in 1:3){
	hist(scale(nda17[,ind_dv[p]]), xlab="", freq=FALSE,main=names(nda17)[ind_dv][p])
}
dev.off()

###################################
###################################
## Subset variables for analyses ##
###################################
###################################

data = nda17[,c(1,ind_nest,ind_demog,which(names(nda17)=="rel_relationship"),which(names(nda17)=="rel_group_id"),ind_np,ind_dv)]
names(data);dim(data)
data$src_subject_id = as.character(data$src_subject_id)
names(data)[names(data)=="src_subject_id"] = "pid"
data$site_num = as.numeric(substr(data$Site,5,6))
data$fam_num = 0
ind=0
for(i in sort(unique(data$rel_family_id))){
	ind = ind+1
	data$fam_num[data$rel_family_id==i & !is.na(data$rel_family_id)] = ind
}
data = data[order(data$site_num,data$fam_num,data$rel_group_id),]

####################
####################
## Fit usual PCA  ## 
####################
####################

## PCA on unimputed data 
ind_Y = c(11:19); names(data)[ind_Y]
Y = as.matrix(scale(data[complete.cases(data[,c(ind_Y)]),ind_Y]))
ev = eigen(cor(Y))
ap = parallel(subject=nrow(Y),var=ncol(Y),rep=100,cent=.05)
nS = nScree(x=ev$values,aparallel=ap$eigen$qevpea)
plotnScree(nS)
ncomp = 3
#y.pca = psych::principal(Y, rotate="promax", nfactors=ncomp, scores=TRUE)
#y.pca$loadings
y.pca = psych::principal(Y, rotate="varimax", nfactors=ncomp, scores=TRUE)
y.pca$loadings

#########################################
#########################################
## Fit 9-variable PCA Model using Stan ##
#########################################
#########################################

###############
## Prep data ##
###############

ind_Y = 11:19; names(data)[ind_Y]
ind_demog = c(4:9); names(data)[ind_demog]
ind_nest = c(2,3); names(data)[ind_nest]

data1 = data[complete.cases(data[,c(ind_Y)]),]; names(data1); dim(data1)
site_num = rep(NA,length(data1$site_num))
for(s in 1:length(unique(data1$site_num))){
	site_num[data1$site_num == unique(data1$site_num)[s]] = s
}
data1$site_num = site_num; rm(site_num)
data1$id_fam = 0
data1$fam_size = 0
ind=0
for(s in 1:length(unique(data1$site_num))){
	data_s = data1[data1$site_num == s, ]
	for(f in 1:length(unique(data_s$rel_family_id))){
		data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
			sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
		if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
			ind=ind+1	
			data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
		}	
	}
	data1[data1$site_num == s, ] = data_s
}
data1 = data1[order(data1$site_num,data1$id_fam),]

Site = data1$site_num
Fam = data1$id_fam
ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
	ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
	nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)
Y = (as.matrix(scale(data1[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

####################
## Run stan model ##
####################

Nsamples = 1000
Nchains = 3

model_file = "stan_code/bppca.stan"
smod = stan_model(model_file)

D_max = 5
sa.list = list()
log_lik.list = list()
looic.list = list()
for(d in 1:D_max){
	print(d)
	pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
	set.seed(314)
	sa.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")
	#log_lik.list[[d]] <- extract_log_lik(sa.list[[d]])
	log_lik.list[[d]] <- extract(sa.list[[d]],"log_lik_marg")[[1]]
	looic.list[[d]] = loo(log_lik.list[[d]])
	save(sa.list,log_lik.list,looic.list,file="results/bppca_results.RData")
	print("###############################")
}	

#################################
## Model selection using LOOIC ##
#################################

load("results/bppca_results.RData")

looic.obj = compare(looic.list[[1]],looic.list[[2]],looic.list[[3]],looic.list[[4]],looic.list[[5]])
print(looic.obj)

d=3

print(sa.list[[d]], pars=c("Q"), probs=c(.025,.5,.975))
print(sa.list[[d]], pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"), probs=c(.025,.5,.975))

## Create SM Figure 3
jpeg("figures_and_tables/fig_sm3.jpeg")
traceplot(sa.list[[d]], pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"))
dev.off()

#######################################
## Extract parameters from the model ##
#######################################

sa = sa.list[[d]]
Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

par(mfrow=c(2,2))
hist(sigma2_a,main="Var(a)",xlim=c(0,.12),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_a,.5),3)))
hist(sigma2_b,main="Var(b)",xlim=c(.35,.65),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_b,.5),3)))
hist(sigma_c^2,main="Var(c)",xlim=c(0,.012),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_c^2,.5),3)))
hist(sigma_d^2,main="Var(d)",xlim=c(.06,.13),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_d^2,.5),3)))

vc_tab = array(0, dim=c(4,3))
vc_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
R_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
S_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
W_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
Theta_old = array(0, dim=c(d,N,Nchains*Nsamples/2))
Theta_new = array(0, dim=c(d,N,Nchains*Nsamples/2))
Lambda_old = array(0, dim=c(P,d,Nchains*Nsamples/2))
Lambda_new = array(0, dim=c(P,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind+ 1
		Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
		Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
		Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
		W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
		R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
		S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
		Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
	}
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
	for(j in 1:(Nsamples*Nchains/2)){
		if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
			ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
			ll_mx[k] = max(abs(Lambda_new[,k,j]))
	}
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		for(k in 1:d){
			Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
		}
		Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
			Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
	}
}

## Varimax & Promax Rotations
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
Lambda_pmx = Lambda_new
Theta_pmx = Theta_new
Rot_pmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		tmp = varimax(Lambda_new[,,ind])
		Rot_vmx[,,ind] = tmp$rot
		Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
		Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
		tmp = psych::Promax(Lambda_new[,,ind])
		Rot_pmx[,,ind] = tmp$rot
		Lambda_pmx[,,ind] = Lambda_new[,,ind]%*%Rot_pmx[,,ind]
		Theta_pmx[,,ind] = t(Rot_pmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
	}
}
## Permute Varimax factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
		Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
		Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
		}
	}
}

## Permute Promax factors
tr=principal(Y,nfactors=d,rotate="promax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_pmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_pmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_pmx[,,ind] = Lambda_pmx[,ll,ind]
		Theta_pmx[,,ind] = Theta_pmx[ll,,ind]
		Rot_pmx[,,ind] = Rot_pmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_pmx[d1,,ind] = Theta_pmx[d1,,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
			Rot_pmx[,d1,ind] = Rot_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
			Lambda_pmx[,d1,ind] = Lambda_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
		}
	}
}
## Variance explained by retained factors
var_new = array(0, dim=c(d,Nchains*Nsamples/2))
var_vmx = array(0, dim=c(d,Nchains*Nsamples/2))
var_pmx = array(0, dim=c(d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
		var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
		var_pmx[,ind] = diag(t(Lambda_pmx[,,ind])%*%(Lambda_pmx[,,ind]))
	}
}
## Communalities & Uniquenesses
comm_uniq_new = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_vmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_pmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		comm_uniq_new[,,ind] = cbind(diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])),1-diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])))
		comm_uniq_vmx[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
		comm_uniq_pmx[,,ind] = cbind(diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])),1-diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])))
	}
}
## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
	for(p2 in 1:P){
		Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
	}
}
comm_uniq_new_m = array(0,dim=c(P,2,3))
comm_uniq_vmx_m = array(0,dim=c(P,2,3))
comm_uniq_pmx_m = array(0,dim=c(P,2,3))
Lambda_new_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
Lambda_pmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
	comm_uniq_new_m[p,1,] = quantile(comm_uniq_new[p,1,],probs=c(.025,.5,.975))
	comm_uniq_vmx_m[p,1,] = quantile(comm_uniq_vmx[p,1,],probs=c(.025,.5,.975))
	comm_uniq_pmx_m[p,1,] = quantile(comm_uniq_pmx[p,1,],probs=c(.025,.5,.975))
	comm_uniq_new_m[p,2,] = quantile(comm_uniq_new[p,2,],probs=c(.025,.5,.975))
	comm_uniq_vmx_m[p,2,] = quantile(comm_uniq_vmx[p,2,],probs=c(.025,.5,.975))
	comm_uniq_pmx_m[p,2,] = quantile(comm_uniq_pmx[p,2,],probs=c(.025,.5,.975))
	for(k in 1:d){
		Lambda_new_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
		Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
		Lambda_pmx_m[p,k,] = quantile(Lambda_pmx[p,k,],probs=c(.025,.5,.975))
	}
}
Theta_new_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
Theta_pmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
	for(i in 1:N){
		Theta_new_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
		Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
		Theta_pmx_m[k,i,] = quantile(Theta_pmx[k,i,],probs=c(.025,.5,.975))
	}
}		
var_new_m = array(0,dim=c(d,3))
var_vmx_m = array(0,dim=c(d,3))
var_pmx_m = array(0,dim=c(d,3))
for(k in 1:d){
	var_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
	var_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
	var_pmx_m[k,] = quantile(var_pmx[k,],c(.025,.5,.975))
}
var_expl_new = apply(var_new_m,2,sum)/P
var_expl_vmx = apply(var_vmx_m,2,sum)/P
var_expl_pmx = apply(var_pmx_m,2,sum)/P

## Create SM Table 1A
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_new_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_new_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_new_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm1a = knitr::kable(tabAsStringMatrix)

## Create SM Table 1B
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_pmx_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_pmx_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_pmx_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm1b = knitr::kable(tabAsStringMatrix)

## Create SM Table 3
comm_uniq.tab = as.data.frame(array(0,dim=c(P+1,6)))
names(comm_uniq.tab) = c("","Communality","","","Uniqueness","")
rownames(comm_uniq.tab) = c("Quantiles:",names(data)[ind_Y])
comm_uniq.tab[1,] = rep(c(".025","0.50",".975"),2)
comm_uniq.tab[2:(P+1),c(1,4)]=sprintf("%.3f", round(comm_uniq_vmx_m[,,1],3))
comm_uniq.tab[2:(P+1),c(2,5)]=sprintf("%.3f", round(comm_uniq_vmx_m[,,2],3))
comm_uniq.tab[2:(P+1),c(3,6)]=sprintf("%.3f", round(comm_uniq_vmx_m[,,3],3))
tabAsStringMatrix = print(comm_uniq.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm3 = knitr::kable(tabAsStringMatrix)

## Create Table 3
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_vmx_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_vmx_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_vmx_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab3 = knitr::kable(tabAsStringMatrix)

## Add PC scores into dataset
data1$pc1 = Theta_vmx_m[1,,2]
data1$pc2 = Theta_vmx_m[2,,2]
data1$pc3 = Theta_vmx_m[3,,2]

##########################
##########################
## Impute Missing Data  ## 
##########################
##########################

data = merge(data,data1[,c(1,which(names(data1)=="pc1"),which(names(data1)=="pc2"),which(names(data1)=="pc3"))],by="pid",all.x=TRUE)

# Number of multiple imputed datasets & maximum number of iterations 
n.imp = 5
n.iter = 5

var.ls <- c("RaceEthnicity","HighestParentalEducation","HouseholdMaritalStatus","HouseholdIncome","PicVocab","Flanker",
				"List","CardSort","Pattern","Picture","Reading","RAVLT","LMT","pc1","pc2","pc3")
dat0 = data.table(data)
dat0 <- dat0[, var.ls, with = FALSE ]
Hmisc::describe(dat0)

# Initialize model
ini = mice( dat0, m = 1, maxit = 0 )
meth = ini$meth
meth["RaceEthnicity"] = "polyreg"
meth["HighestParentalEducation"] = "polyreg"
meth["HouseholdMaritalStatus"] = "polyreg"
meth["HouseholdIncome"] = "polyreg"
meth["PicVocab"] = "norm.predict"
meth["Flanker"] = "norm.predict"
meth["List"] = "norm.predict"
meth["CardSort"] = "norm.predict"
meth["Pattern"] = "norm.predict"
meth["Picture"] = "norm.predict"
meth["Reading"] = "norm.predict"
meth["RAVLT"] = "norm.predict"
meth["LMT"] = "norm.predict"
meth["pc1"] = "norm.predict"
meth["pc2"] = "norm.predict"
meth["pc3"] = "norm.predict"
pred = ini$pred

set.seed(314)
post = mice( dat0, meth = meth, pred = pred, seed = 111,
              m = 1, maxit = 0)$post
dat.imp = mice( dat0, meth = meth, pred = pred, post = post,
                 seed = 1111, m = n.imp, maxit = n.iter)
rm(dat0)
# Convert imputed data to long format
dat.mi <- complete(dat.imp, action = "long", include = TRUE)
rm(dat.imp)

names(dat.mi)[1:2] <- c("imp", "id")
dat.mi$pid = rep(data$pid,n.imp+1)

# Adding variables
dat = data[,c(1,2,3,4,5,20:22)]
dat.mi <- merge( dat.mi, dat, by = "pid", all.x = TRUE )
dat.mi$id = NULL
names(dat.mi)

######################################################
######################################################
## Association of factor scores with CBCL variables ##
######################################################
######################################################

d_rsq = array(NA, dim = c(n.imp, 3))
gamm.results = list()
for(j in 1:n.imp){
	dat.j = dat.mi[dat.mi$imp==j,]
	ind_y = c(which(names(dat.j)=="Externalizing"),which(names(dat.j)=="Internalizing"),which(names(dat.j)=="Stress")); names(dat.j)[ind_y]
	dat.j[,ind_y] = dat.j[,ind_y]+.1
	## Run GAMMs
	p.values = array(NA, dim=c(1,length(ind_y)))
	t.values = array(NA, dim=c(1,length(ind_y)))
	results = list()
	for(i in 1:length(ind_y)){
		show(c(j,i))
		results[[i]] = list()
		form0 = paste(names(dat.j)[ind_y[i]],"~ ")
		for(k in 1:length(ind_demog)){
			form0 = paste(form0,"+",names(data1)[ind_demog[k]])
		}
		form0 = formula(form0)	
		form = paste(names(dat.j)[ind_y[i]],"~ pc1+pc2+pc3")
		for(k in 1:length(ind_demog)){
			form = paste(form,"+",names(data1)[ind_demog[k]])
		}
		form = formula(form)	
		ran = paste("~(1|",names(data1)[ind_nest[1]])
		for(k in 2:length(ind_nest)){
			ran = paste(ran,"/",names(data1)[ind_nest[k]])
		}	
		ran = paste(ran,")")
		ran = formula(ran)
		gamm0 = gamm4(formula = form0, random = ran, data = dat.j, family = Gamma(link = "log"))
		gamm1 = gamm4(formula = form, random = ran, data = dat.j, family = Gamma(link = "log"))
		results[[i]][[1]] = gamm0
		results[[i]][[2]] = gamm1
		p.values[1,i] = summary(results[[i]][[2]]$gam)$p.pv[2]
		t.values[1,i] = summary(results[[i]][[2]]$gam)$p.t[2]
		print(summary(results[[i]][[1]]$gam))
		print(summary(results[[i]][[2]]$gam))
		print(anova(results[[i]][[2]]$gam))
		print("######################################")
	}
	d_rsq[j,] = rep(0,length(ind_y))
	for(p in 1:length(ind_y)){
		d_rsq[j,p] = summary(results[[p]][[2]]$gam)$r.sq - summary(results[[p]][[1]]$gam)$r.sq
	}
	gamm.results[[j]] = results
}	
d_rsq = apply(d_rsq,2,mean)

## Create Table 4
cbcl.gamm.tab = as.data.frame(array(0,dim=c(15+d,length(ind_y)*3)))
rownames(cbcl.gamm.tab) = c("Variable",names(summary(results[[1]][[2]]$gam)$p.coeff))
names(cbcl.gamm.tab) = c(names(data1)[ind_y][1],names(data1)[ind_y][1],names(data1)[ind_y][1],
							names(data1)[ind_y][2],names(data1)[ind_y][2],names(data1)[ind_y][2],
							names(data1)[ind_y][3],names(data1)[ind_y][3],names(data1)[ind_y][3])
cbcl.gamm.tab[1,] = rep(c("coef   ","se    ","p_value"),length(ind_y))
ind = 0
for(k in 1:length(ind_y)){
	ind = ind + 1
	coeff = summary(gamm.results[[1]][[k]][[2]]$gam)$p.table[,1]
	for(b in 2:n.imp){
		coeff = 	coeff + summary(gamm.results[[b]][[k]][[2]]$gam)$p.table[,1]
	}
	coeff = coeff/n.imp	
	cbcl.gamm.tab[2:(15+d),ind] = coeff
	ind = ind + 1
	se = summary(gamm.results[[1]][[k]][[2]]$gam)$p.table[,2]
	for(b in 2:n.imp){
		se = se + summary(gamm.results[[b]][[k]][[2]]$gam)$p.table[,2]
	}
	se = se/n.imp
	imp.err = 0*se
	for(b in 1:n.imp){
		imp.err = imp.err + (summary(gamm.results[[b]][[k]][[2]]$gam)$p.table[,1]-coeff)^2
	}
	imp.err = imp.err/n.imp	
	se = se + imp.err
	cbcl.gamm.tab[2:(15+d),ind] = round(se,2)
	ind = ind + 1
	cbcl.gamm.tab[2:(15+d),ind] = 2*(1-pnorm(abs(coeff),0,se))
}
cbcl.gamm.tab
tabAsStringMatrix = print(cbcl.gamm.tab, printToggle = FALSE, noSpaces = TRUE)
tab4 = knitr::kable(tabAsStringMatrix)

#############################
## K-fold Cross Validation ##
#############################

K = 10

set.seed(314)
cv = cvFolds(n=dim(data)[1], K = K, R = 1, type = c("random"))

gamm.cv.results = list()
for(k in 1:K){
	dat.k = data[!cv$which==k,]
	dat.test = data[cv$which==k,]
	ind_y = c(which(names(dat.k)=="Externalizing"),which(names(dat.k)=="Internalizing"),which(names(dat.k)=="Stress")); names(dat.k)[ind_y]
	dat.k[,ind_y] = dat.k[,ind_y]+.1
	## Run GAMMs
	results = list()
	for(i in 1:length(ind_y)){
		show(c(k,i))
		results[[i]] = list()
		form0 = paste(names(dat.k)[ind_y[i]],"~ ")
		for(kk in 1:length(ind_demog)){
			form0 = paste(form0,"+",names(data1)[ind_demog[kk]])
		}
		form0 = formula(form0)	
		form = paste(names(dat.k)[ind_y[i]],"~ pc1+pc2+pc3")
		for(kk in 1:length(ind_demog)){
			form = paste(form,"+",names(data1)[ind_demog[kk]])
		}
		form = formula(form)	
		ran = paste("~(1|",names(data1)[ind_nest[1]])
		for(kk in 2:length(ind_nest)){
			ran = paste(ran,"/",names(data1)[ind_nest[kk]])
		}	
		ran = paste(ran,")")
		ran = formula(ran)
		gamm0 = gamm4(formula = form0, random = ran, data = dat.k, family = Gamma(link = "log"))
		gamm1 = gamm4(formula = form, random = ran, data = dat.k, family = Gamma(link = "log"))
		results[[i]][[1]] = gamm0
		results[[i]][[2]] = gamm1
		print(summary(results[[i]][[1]]$gam))
		print(summary(results[[i]][[2]]$gam))
		print("######################################")
	}
	gamm.cv.results[[k]] = results
}	

d_rsq.cv = array(0, dim = c(K, 3))
for(k in 1:K){
	results = gamm.cv.results[[k]]
	dat.k = data[!cv$which==k,]
	dat.test = data[cv$which==k,]	
	for(p in 1:length(ind_y)){
		d_rsq.cv[k,p] = cor.test(predict(results[[p]][[2]]$gam, newdata = dat.test),dat.test[,ind_y[p]],method = "spearman")$estimate^2 -
					 	cor.test(predict(results[[p]][[1]]$gam, newdata = dat.test),dat.test[,ind_y[p]],method = "spearman")$estimate^2
	}
}
d_rsq.cv = apply(d_rsq.cv,2,mean)

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

##########################################
##########################################
## Fit NIH Toolbox PCA Model using Stan ##
##########################################
##########################################

###############
## Prep data ##
###############

ind_Y = c(11:17); names(data)[ind_Y]
ind_demog = c(4:9); names(data)[ind_demog]
ind_nest = c(2,3); names(data)[ind_nest]

data2 = data[complete.cases(data[,c(ind_Y,ind_nest)]),]; names(data2); dim(data2)
site_num = rep(NA,length(data2$site_num))
for(s in 1:length(unique(data2$site_num))){
	site_num[data2$site_num == unique(data2$site_num)[s]] = s
}
data2$site_num = site_num; rm(site_num)
data2$id_fam = 0
data2$fam_size = 0
ind=0
for(s in 1:length(unique(data2$site_num))){
	data_s = data2[data2$site_num == s, ]
	for(f in 1:length(unique(data_s$rel_family_id))){
		data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
			sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
		if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
			ind=ind+1	
			data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
		}	
	}
	data2[data2$site_num == s, ] = data_s
}
data2 = data2[order(data2$site_num,data2$id_fam),]

Site = data2$site_num
Fam = data2$id_fam

ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
	ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
	nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)

Y = (as.matrix(scale(data2[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

####################
## Run stan model ##
####################

Nsamples = 1000
Nchains = 3

model_file = "stan_code/bppca.stan"
smod = stan_model(model_file)

D_max = 4
sa.nihtb.list = list()
log_lik.nihtb.list = list()
looic.nihtb.list = list()
set.seed(314)
for(d in 1:D_max){
	print(d)
	pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
	set.seed(314)
	sa.nihtb.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")  
	log_lik.nihtb.list[[d]] <- extract(sa.nihtb.list[[d]],"log_lik_marg")[[1]]
	looic.nihtb.list[[d]] = loo(log_lik.nihtb.list[[d]])
	save(sa.nihtb.list,log_lik.nihtb.list,looic.nihtb.list,file="results/bppca_nihtb_results.RData")
	print("###############################")
}

#################################
## Model selection using LOOIC ##
#################################

load("results/bppca_nihtb_results.Rdata")

looic.nihtb.obj = compare(looic.nihtb.list[[1]],looic.nihtb.list[[2]],looic.nihtb.list[[3]],looic.nihtb.list[[4]])
print(looic.nihtb.obj)

d=3

## Create SM Figure 4
pdf("figures_and_tables/fig_sm4.pdf")
traceplot(sa.nihtb.list[[d]], pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"))
dev.off()


#######################################
## Extract parameters from the model ##
#######################################

sa = sa.nihtb.list[[d]]
Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

par(mfrow=c(2,2))
hist(sigma2_a,main="Var(a)",xlim=c(0,.12),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_a,.5),3)))
hist(sigma2_b,main="Var(b)",xlim=c(.35,.65),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_b,.5),3)))
hist(sigma_c^2,main="Var(c)",xlim=c(0,.012),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_c^2,.5),3)))
hist(sigma_d^2,main="Var(d)",xlim=c(.06,.13),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_d^2,.5),3)))

vc_nihtb_tab = array(0, dim=c(4,3))
vc_nihtb_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_nihtb_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_nihtb_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_nihtb_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
R_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
S_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
W_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
Theta_old = array(0, dim=c(d,N,Nchains*Nsamples/2))
Theta_new = array(0, dim=c(d,N,Nchains*Nsamples/2))
Lambda_old = array(0, dim=c(P,d,Nchains*Nsamples/2))
Lambda_new = array(0, dim=c(P,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind+ 1
		Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
		Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
		Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
		W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
		R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
		S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
		Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
	}
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
	for(j in 1:(Nsamples*Nchains/2)){
		if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
			ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
			ll_mx[k] = max(abs(Lambda_new[,k,j]))
	}
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		for(k in 1:d){
			Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
		}
		Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
			Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
	}
}

## Varimax & Promax Rotations
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
Lambda_pmx = Lambda_new
Theta_pmx = Theta_new
Rot_pmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		tmp = varimax(Lambda_new[,,ind])
		Rot_vmx[,,ind] = tmp$rot
		Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
		Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
		tmp = psych::Promax(Lambda_new[,,ind])
		Rot_pmx[,,ind] = tmp$rot
		Lambda_pmx[,,ind] = Lambda_new[,,ind]%*%Rot_pmx[,,ind]
		Theta_pmx[,,ind] = t(Rot_pmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
	}
}

## Permute Varimax factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
tr = tr[,c(2,1,3)]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
		Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
		Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
		}
	}
}

## Permute Promax factors
tr=principal(Y,nfactors=d,rotate="promax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_pmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_pmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_pmx[,,ind] = Lambda_pmx[,ll,ind]
		Theta_pmx[,,ind] = Theta_pmx[ll,,ind]
		Rot_pmx[,,ind] = Rot_pmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_pmx[d1,,ind] = Theta_pmx[d1,,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
			Rot_pmx[,d1,ind] = Rot_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
			Lambda_pmx[,d1,ind] = Lambda_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
		}
	}
}

## Variance explained by retained factors
var_new = array(0, dim=c(d,Nchains*Nsamples/2))
var_vmx = array(0, dim=c(d,Nchains*Nsamples/2))
var_pmx = array(0, dim=c(d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
		var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
		var_pmx[,ind] = diag(t(Lambda_pmx[,,ind])%*%(Lambda_pmx[,,ind]))
	}
}
## Communalities & Uniquenesses
comm_uniq_new = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_vmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_pmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		comm_uniq_new[,,ind] = cbind(diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])),1-diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])))
		comm_uniq_vmx[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
		comm_uniq_pmx[,,ind] = cbind(diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])),1-diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])))
	}
}
## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
	for(p2 in 1:P){
		Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
	}
}
comm_uniq_new_m = array(0,dim=c(P,2,3))
comm_uniq_vmx_m = array(0,dim=c(P,2,3))
comm_uniq_pmx_m = array(0,dim=c(P,2,3))
Lambda_new_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
Lambda_pmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
	comm_uniq_new_m[p,1,] = quantile(comm_uniq_new[p,1,],probs=c(.025,.5,.975))
	comm_uniq_vmx_m[p,1,] = quantile(comm_uniq_vmx[p,1,],probs=c(.025,.5,.975))
	comm_uniq_pmx_m[p,1,] = quantile(comm_uniq_pmx[p,1,],probs=c(.025,.5,.975))
	comm_uniq_new_m[p,2,] = quantile(comm_uniq_new[p,2,],probs=c(.025,.5,.975))
	comm_uniq_vmx_m[p,2,] = quantile(comm_uniq_vmx[p,2,],probs=c(.025,.5,.975))
	comm_uniq_pmx_m[p,2,] = quantile(comm_uniq_pmx[p,2,],probs=c(.025,.5,.975))
	for(k in 1:d){
		Lambda_new_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
		Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
		Lambda_pmx_m[p,k,] = quantile(Lambda_pmx[p,k,],probs=c(.025,.5,.975))
	}
}
Theta_new_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
Theta_pmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
	for(i in 1:N){
		Theta_new_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
		Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
		Theta_pmx_m[k,i,] = quantile(Theta_pmx[k,i,],probs=c(.025,.5,.975))
	}
}		
var_nihtb_new_m = array(0,dim=c(d,3))
var_nihtb_vmx_m = array(0,dim=c(d,3))
var_nihtb_pmx_m = array(0,dim=c(d,3))
for(k in 1:d){
	var_nihtb_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
	var_nihtb_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
	var_nihtb_pmx_m[k,] = quantile(var_pmx[k,],c(.025,.5,.975))
}
var_nihtb_expl_new = apply(var_new_m,2,sum)/P
var_nihtb_expl_vmx = apply(var_vmx_m,2,sum)/P
var_nihtb_expl_pmx = apply(var_pmx_m,2,sum)/P

## Create SM Table 2A
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_new_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_new_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_new_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm2a = knitr::kable(tabAsStringMatrix)

## Create SM Table 2B
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_pmx_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_pmx_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_pmx_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm2b = knitr::kable(tabAsStringMatrix)

## Create SM Table 4
comm_uniq.nihtb.tab = as.data.frame(array(0,dim=c(P+1,6)))
names(comm_uniq.nihtb.tab) = c("","Communality","","","Uniqueness","")
rownames(comm_uniq.nihtb.tab) = c("Quantiles:",names(data)[ind_Y])
comm_uniq.nihtb.tab[1,] = rep(c(".025","0.50",".975"),2)
comm_uniq.nihtb.tab[2:(P+1),c(1,4)]=sprintf("%.3f", round(comm_uniq_new_m[,,1],3))
comm_uniq.nihtb.tab[2:(P+1),c(2,5)]=sprintf("%.3f", round(comm_uniq_new_m[,,2],3))
comm_uniq.nihtb.tab[2:(P+1),c(3,6)]=sprintf("%.3f", round(comm_uniq_new_m[,,3],3))
tabAsStringMatrix = print(comm_uniq.nihtb.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm4 = knitr::kable(tabAsStringMatrix)

## Create Table 5
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_vmx_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_vmx_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_vmx_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab5 = knitr::kable(tabAsStringMatrix)

data2$pc1 = Theta_vmx_m[1,,2]
data2$pc2 = Theta_vmx_m[2,,2]
data2$pc3 = Theta_vmx_m[3,,2]

##########################
##########################
## Impute Missing Data  ## 
##########################
##########################

data$pc1=NULL
data$pc2=NULL
data$pc3=NULL

data = merge(data,data2[,c(1,which(names(data2)=="pc1"),which(names(data2)=="pc2"),which(names(data2)=="pc3"))],by="pid",all.x=TRUE)

# Number of multiple imputed datasets & maximum number of iterations 
n.imp = 5
n.iter = 5

var.ls <- c("RaceEthnicity","HighestParentalEducation","HouseholdMaritalStatus","HouseholdIncome","PicVocab","Flanker",
				"List","CardSort","Pattern","Picture","Reading","pc1","pc2","pc3")
dat0 = data.table(data)
dat0 <- dat0[, var.ls, with = FALSE ]
Hmisc::describe(dat0)

# Initialize model
ini = mice( dat0, m = 1, maxit = 0 )
meth = ini$meth
meth["RaceEthnicity"] = "polyreg"
meth["HighestParentalEducation"] = "polyreg"
meth["HouseholdMaritalStatus"] = "polyreg"
meth["HouseholdIncome"] = "polyreg"
meth["PicVocab"] = "norm.predict"
meth["Flanker"] = "norm.predict"
meth["List"] = "norm.predict"
meth["CardSort"] = "norm.predict"
meth["Pattern"] = "norm.predict"
meth["Picture"] = "norm.predict"
meth["Reading"] = "norm.predict"
meth["pc1"] = "norm.predict"
meth["pc2"] = "norm.predict"
meth["pc3"] = "norm.predict"
pred = ini$pred

set.seed(314)
post = mice( dat0, meth = meth, pred = pred, seed = 111,
              m = 1, maxit = 0)$post
dat.imp = mice( dat0, meth = meth, pred = pred, post = post,
                 seed = 1111, m = n.imp, maxit = n.iter)
rm(dat0)
# Convert imputed data to long format
dat.mi <- complete(dat.imp, action = "long", include = TRUE)
rm(dat.imp)

names(dat.mi)[1:2] <- c("imp", "id")
dat.mi$pid = rep(data$pid,n.imp+1)

# Adding variables
dat = data[,c(1,2,3,4,5,20:22)]; names(dat)
dat.mi <- merge( dat.mi, dat, by = "pid", all.x = TRUE )
dat.mi$id = NULL
names(dat.mi)

######################################################
######################################################
## Association of factor scores with CBCL variables ##
######################################################
######################################################

d_rsq_nihtbx = array(0, dim = c(n.imp, 3))
gamm.results.nihtbx = list()
for(j in 1:n.imp){
	dat.j = dat.mi[dat.mi$imp==j,]
	ind_y = c(which(names(dat.j)=="Externalizing"),which(names(dat.j)=="Internalizing"),which(names(dat.j)=="Stress")); names(dat.j)[ind_y]
	dat.j[,ind_y] = dat.j[,ind_y]+.1
	## Run GAMMs
	p.values = array(NA, dim=c(1,length(ind_y)))
	t.values = array(NA, dim=c(1,length(ind_y)))
	results = list()
	for(i in 1:length(ind_y)){
		show(c(j,i))
		results[[i]] = list()
		form0 = paste(names(dat.j)[ind_y[i]],"~ ")
		for(k in 1:length(ind_demog)){
			form0 = paste(form0,"+",names(data1)[ind_demog[k]])
		}
		form0 = formula(form0)	
		form = paste(names(dat.j)[ind_y[i]],"~ pc1+pc2+pc3")
		for(k in 1:length(ind_demog)){
			form = paste(form,"+",names(data1)[ind_demog[k]])
		}
		form = formula(form)	
		ran = paste("~(1|",names(data1)[ind_nest[1]])
		for(k in 2:length(ind_nest)){
			ran = paste(ran,"/",names(data1)[ind_nest[k]])
		}	
		ran = paste(ran,")")
		ran = formula(ran)
		gamm0 = gamm4(formula = form0, random = ran, data = dat.j, family = Gamma(link = "log"))
		gamm1 = gamm4(formula = form, random = ran, data = dat.j, family = Gamma(link = "log"))
		results[[i]][[1]] = gamm0
		results[[i]][[2]] = gamm1
		p.values[1,i] = summary(results[[i]][[2]]$gam)$p.pv[2]
		t.values[1,i] = summary(results[[i]][[2]]$gam)$p.t[2]
		print(summary(results[[i]][[1]]$gam))
		print(summary(results[[i]][[2]]$gam))
		print(anova(results[[i]][[2]]$gam))
		print("######################################")
	}
	for(p in 1:length(ind_y)){
		d_rsq_nihtbx[j,p] = summary(results[[p]][[2]]$gam)$r.sq - summary(results[[p]][[1]]$gam)$r.sq
	}
	gamm.results.nihtbx[[j]] = results
}	
d_rsq_nihtbx = apply(d_rsq_nihtbx,2,mean)


## Create Table 6
cbcl.gamm.tab = as.data.frame(array(0,dim=c(15+d,length(ind_y)*3)))
rownames(cbcl.gamm.tab) = c("Variable",names(summary(results[[1]][[2]]$gam)$p.coeff))
names(cbcl.gamm.tab) = c(names(data1)[ind_y][1],names(data1)[ind_y][1],names(data1)[ind_y][1],
							names(data1)[ind_y][2],names(data1)[ind_y][2],names(data1)[ind_y][2],
							names(data1)[ind_y][3],names(data1)[ind_y][3],names(data1)[ind_y][3])
cbcl.gamm.tab[1,] = rep(c("coef   ","se    ","p_value"),length(ind_y))
ind = 0
for(k in 1:length(ind_y)){
	ind = ind + 1
	coeff = summary(gamm.results.nihtbx[[1]][[k]][[2]]$gam)$p.table[,1]
	for(b in 2:n.imp){
		coeff = 	coeff + summary(gamm.results.nihtbx[[b]][[k]][[2]]$gam)$p.table[,1]
	}
	coeff = coeff/n.imp	
	cbcl.gamm.tab[2:(15+d),ind] = coeff
	ind = ind + 1
	se = summary(gamm.results.nihtbx[[1]][[k]][[2]]$gam)$p.table[,2]
	for(b in 2:n.imp){
		se = se + summary(gamm.results.nihtbx[[b]][[k]][[2]]$gam)$p.table[,2]
	}
	se = se/n.imp
	imp.err = 0*se
	for(b in 1:n.imp){
		imp.err = imp.err + (summary(gamm.results.nihtbx[[b]][[k]][[2]]$gam)$p.table[,1]-coeff)^2
	}
	imp.err = imp.err/n.imp	
	se = se + imp.err
	cbcl.gamm.tab[2:(15+d),ind] = round(se,2)
	ind = ind + 1
	cbcl.gamm.tab[2:(15+d),ind] = 2*(1-pnorm(abs(coeff),0,se))
}
cbcl.gamm.tab
tabAsStringMatrix = print(cbcl.gamm.tab, printToggle = FALSE, noSpaces = TRUE)
tab6 = knitr::kable(tabAsStringMatrix)


#############################
## K-fold Cross Validation ##
#############################

K = 10

set.seed(314)
cv = cvFolds(n=dim(data)[1], K = K, R = 1, type = c("random"))

gamm.cv.nihtbx.results = list()
for(k in 1:K){
	dat.k = data[!cv$which==k,]
	dat.test = data[cv$which==k,]
	ind_y = c(which(names(dat.k)=="Externalizing"),which(names(dat.k)=="Internalizing"),which(names(dat.k)=="Stress")); names(dat.k)[ind_y]
	dat.k[,ind_y] = dat.k[,ind_y]+.1
	## Run GAMMs
	results = list()
	for(i in 1:length(ind_y)){
		show(c(k,i))
		results[[i]] = list()
		form0 = paste(names(dat.k)[ind_y[i]],"~ ")
		for(kk in 1:length(ind_demog)){
			form0 = paste(form0,"+",names(data1)[ind_demog[kk]])
		}
		form0 = formula(form0)	
		form = paste(names(dat.k)[ind_y[i]],"~ pc1+pc2+pc3")
		for(kk in 1:length(ind_demog)){
			form = paste(form,"+",names(data1)[ind_demog[kk]])
		}
		form = formula(form)	
		ran = paste("~(1|",names(data1)[ind_nest[1]])
		for(kk in 2:length(ind_nest)){
			ran = paste(ran,"/",names(data1)[ind_nest[kk]])
		}	
		ran = paste(ran,")")
		ran = formula(ran)
		gamm0 = gamm4(formula = form0, random = ran, data = dat.k, family = Gamma(link = "log"))
		gamm1 = gamm4(formula = form, random = ran, data = dat.k, family = Gamma(link = "log"))
		results[[i]][[1]] = gamm0
		results[[i]][[2]] = gamm1
		print(summary(results[[i]][[1]]$gam))
		print(summary(results[[i]][[2]]$gam))
		print("######################################")
	}
	gamm.cv.nihtbx.results[[k]] = results
}	

d_rsq.cv.nihtbx = array(0, dim = c(K, 3))
for(k in 1:K){
	results = gamm.cv.nihtbx.results[[k]]
	dat.k = data[!cv$which==k,]
	dat.test = data[cv$which==k,]	
	for(p in 1:length(ind_y)){
		d_rsq.cv.nihtbx[k,p] = cor.test(predict(results[[p]][[2]]$gam, newdata = dat.test),dat.test[,ind_y[p]],method = "spearman")$estimate^2 -
					 	cor.test(predict(results[[p]][[1]]$gam, newdata = dat.test),dat.test[,ind_y[p]],method = "spearman")$estimate^2
	}
}
d_rsq.cv.nihtbx = apply(d_rsq.cv.nihtbx,2,mean)


######################
######################
## Create Figure 1  ##
######################
######################

data2$nihtb_pc1 = data2$pc1
data2$nihtb_pc2 = data2$pc2
data2$nihtb_pc3 = data2$pc3
data3  = merge(data1,data2[,c(1,grep("nihtb_pc",names(data2)))], by="pid")
vars =c("Externalizing","Internalizing","Stress","pc1","pc2","pc3","nihtb_pc1","nihtb_pc2","nihtb_pc3")
vars = vars[vars%in% names(data3)]
vars = data3[,vars]
cormat = cor(vars , use="pairwise.complete.obs", method="spearman")
cormat <- rcorr(as.matrix(vars), type="spearman")
melt_cormat = melt(cormat$r)
melt_cormat2 = melt_cormat
melt_cormat2[melt_cormat[,3]==1,3]=NA
for(j in 1:81){
	tmp = melt_cormat2[melt_cormat[,3]	== melt_cormat[j,3],3]
	if(length(tmp)==2) tmp[1] = NA
	melt_cormat2[melt_cormat[,3]	== melt_cormat[j,3],3]	= tmp 
}
melt_cormat2$value = round(melt_cormat2$value, 2)
melt_cormat = melt_cormat2; rm(melt_cormat2)

pdf("figures_and_tables/fig1.pdf")
ggplot(data = melt_cormat, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = "white")+
  xlab("") + ylab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust = 1)) + 
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 2.5) 
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#############################
#############################
## Split-half BPPCA models ##
#############################
#############################

set.seed(314)
t1 = sample(dim(data)[1],dim(data)[1]/2)
data_t1 = data[t1,]
data_t2 = data[-t1,]

####################################
## Prep data for first split half ##
####################################

ind_Y = c(11:19); names(data_t1)[ind_Y]
ind_demog = 4:9; names(data_t1)[ind_demog]
ind_nest = c(2,3); names(data_t1)[ind_nest]

data11 = data_t1[complete.cases(data_t1[,c(ind_Y,ind_nest)]),]; names(data11); dim(data11)
site_num = rep(NA,length(data11$site_num))
for(s in 1:length(unique(data11$site_num))){
	site_num[data11$site_num == unique(data11$site_num)[s]] = s
}
data11$site_num = site_num; rm(site_num)
data11$id_fam = 0
data11$fam_size = 0
ind=0
for(s in 1:length(unique(data11$site_num))){
	data_s = data11[data11$site_num == s, ]
	for(f in 1:length(unique(data_s$rel_family_id))){
		data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
			sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
		if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
			ind=ind+1	
			data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
		}	
	}
	data11[data11$site_num == s, ] = data_s
}
data11 = data11[order(data11$site_num,data11$id_fam),]

Site = data11$site_num
Fam = data11$id_fam
ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
	ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
	nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)
Y = (as.matrix(scale(data11[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

####################
## Run stan model ##
####################

Nsamples = 1000
Nchains = 3

set.seed(314)
model_file = "stan_code/bppca.stan"
smod = stan_model(model_file)

D_max = 3
sa.splithalf1.list = list()
log_lik.splithalf1.list = list()
looic.splithalf1.list = list()
for(d in D_max){
	print(d)
	pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
	set.seed(314)
	sa.splithalf1.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")  
	log_lik.splithalf1.list[[d]] <- extract(sa.splithalf1.list[[d]],"log_lik_marg")[[1]]
	looic.splithalf1.list[[d]] = loo(log_lik.splithalf1.list[[d]])
	save(sa.splithalf1.list,log_lik.splithalf1.list,
		looic.splithalf1.list,file="results/bppca_splithalf1_results.RData")
	print("###############################")
}

#looic.splithalf1.obj = compare(looic.splithalf1.list[[1]],looic.splithalf1.list[[2]],looic.splithalf1.list[[3]],
#	looic.splithalf1.list[[4]],looic.splithalf1.list[[5]],looic.splithalf1.list[[6]],
#	looic.splithalf1.list[[7]],looic.splithalf1.list[[8]])
#print(looic.splithalf1.obj)


###################################
## Extract parameters from model ##
###################################

load("results/bppca_splithalf1_results.RData")

d = 3
sa = sa.splithalf1.list[[d]]

Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

print(sa, pars=c("Q"), probs=c(.025,.5,.975))
print(sa, pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"), probs=c(.025,.5,.975))
traceplot(sa, pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"))

vc_splithalf1_tab = array(0, dim=c(4,3))
vc_splithalf1_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_splithalf1_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_splithalf1_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_splithalf1_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nsamples*Nchains/2))
R_new = array(0, dim=c(P,P,Nsamples*Nchains/2))
S_new = array(0, dim=c(d,d,Nsamples*Nchains/2)) 
W_new = array(0, dim=c(d,d,Nsamples*Nchains/2)) 
Theta_old = array(0, dim=c(d,N,Nsamples*Nchains/2))
Theta_new = array(0, dim=c(d,N,Nsamples*Nchains/2))
Lambda_old = array(0, dim=c(P,d,Nsamples*Nchains/2))
Lambda_new = array(0, dim=c(P,d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind+ 1
		Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
		Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
		Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
		W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
		R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
		S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
		Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
	}
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
	for(j in 1:(Nsamples*Nchains/2)){
		if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
			ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
			ll_mx[k] = max(abs(Lambda_new[,k,j]))
	}
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		for(k in 1:d){
			Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
		}
		Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
			Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
	}
}

## Varimax Rotation
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		tmp = varimax(Lambda_new[,,ind])
		Rot_vmx[,,ind] = tmp$rot
		Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
		Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
	}
}
## Permute factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
		Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
		Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
		}
	}
}

## Variance explained by retained factors
var_new = array(0, dim=c(d,Nsamples*Nchains/2))
var_vmx = array(0, dim=c(d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
		var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
	}
}

## Communalities & Uniquenesses
comm_uniq = array(0, dim=c(P,2,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		comm_uniq[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
	}
}

## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
	for(p2 in 1:P){
		Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
	}
}
comm_uniq_m = array(0,dim=c(P,2,3))
Lambda_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
	comm_uniq_m[p,1,] = quantile(comm_uniq[p,1,],probs=c(.025,.5,.975))
	comm_uniq_m[p,2,] = quantile(comm_uniq[p,2,],probs=c(.025,.5,.975))
	for(k in 1:d){
		Lambda_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
		Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
	}
}
Theta_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
	for(i in 1:N){
		Theta_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
		Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
	}
}		
var_splithalf1_new_m = array(0,dim=c(d,3))
var_splithalf1_vmx_m = array(0,dim=c(d,3))
for(k in 1:d){
	var_splithalf1_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
	var_splithalf1_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
}
var_splithalf1_expl = apply(var_splithalf1_new_m,2,sum)/P
var_splithalf1_vmx_expl = apply(var_splithalf1_vmx_m,2,sum)/P

lambda_splithalf1.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda_splithalf1.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda_splithalf1.tab) = c("",names(data)[ind_Y])
lambda_splithalf1.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda_splithalf1.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_vmx_m[,d1,1],3)
	lambda_splithalf1.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_vmx_m[,d1,2],3)
	lambda_splithalf1.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_vmx_m[,d1,3],3)
}

## Create SM Table 5
tabAsStringMatrix = print(lambda_splithalf1.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm5 = knitr::kable(tabAsStringMatrix)

###############################
## Prep data for second half ##
###############################

ind_Y = c(11:19); names(data_t2)[ind_Y]
ind_demog = 4:9; names(data_t2)[ind_demog]
ind_nest = c(2,3); names(data_t2)[ind_nest]

data12 = data_t2[complete.cases(data_t2[,c(ind_Y,ind_nest)]),]; names(data12); dim(data12)
site_num = rep(NA,length(data12$site_num))
for(s in 1:length(unique(data12$site_num))){
	site_num[data12$site_num == unique(data12$site_num)[s]] = s
}
data12$site_num = site_num; rm(site_num)
data12$id_fam = 0
data12$fam_size = 0
ind=0
for(s in 1:length(unique(data12$site_num))){
	data_s = data12[data12$site_num == s, ]
	for(f in 1:length(unique(data_s$rel_family_id))){
		data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
			sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
		if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
			ind=ind+1	
			data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
		}	
	}
	data12[data12$site_num == s, ] = data_s
}
data12 = data12[order(data12$site_num,data12$id_fam),]

Site = data12$site_num
Fam = data12$id_fam
ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
	ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
	nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)
Y = (as.matrix(scale(data12[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

###################################
## Run stan model on second half ##
###################################

Nsamples = 1000
Nchains = 3

set.seed(314)
model_file = "stan_code/bppca.stan"
smod = stan_model(model_file)

D_max = 3
sa.splithalf2.list = list()
log_lik.splithalf2.list = list()
looic.splithalf2.list = list()
for(d in D_max){
	print(d)
	pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
	set.seed(314)
	sa.splithalf2.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")  
	log_lik.splithalf2.list[[d]] <- extract(sa.splithalf2.list[[d]],"log_lik_marg")[[1]]
	looic.splithalf2.list[[d]] = loo(log_lik.splithalf2.list[[d]])
	save(sa.splithalf2.list,log_lik.splithalf2.list,
		looic.splithalf2.list,file="results/bppca_splithalf2_results.RData")
	print("###############################")
}

#looic.splithalf2.obj = compare(looic.splithalf2.list[[1]],looic.splithalf2.list[[2]],looic.splithalf2.list[[3]],
#	looic.splithalf2.list[[4]],looic.splithalf2.list[[5]],looic.splithalf2.list[[6]],
#	looic.splithalf2.list[[7]],looic.splithalf2.list[[8]])
#print(looic.splithalf2.obj)

###################################
## Extract parameters from model ##
###################################

load("results/bppca_splithalf2_results.RData")

d = 3
sa = sa.splithalf2.list[[d]]

Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

print(sa, pars=c("Q"), probs=c(.025,.5,.975))
print(sa, pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"), probs=c(.025,.5,.975))
traceplot(sa, pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"))

vc_splithalf2_tab = array(0, dim=c(4,3))
vc_splithalf2_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_splithalf2_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_splithalf2_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_splithalf2_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nsamples*Nchains/2))
R_new = array(0, dim=c(P,P,Nsamples*Nchains/2))
S_new = array(0, dim=c(d,d,Nsamples*Nchains/2)) 
W_new = array(0, dim=c(d,d,Nsamples*Nchains/2)) 
Theta_old = array(0, dim=c(d,N,Nsamples*Nchains/2))
Theta_new = array(0, dim=c(d,N,Nsamples*Nchains/2))
Lambda_old = array(0, dim=c(P,d,Nsamples*Nchains/2))
Lambda_new = array(0, dim=c(P,d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind+ 1
		Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
		Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
		Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
		W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
		R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
		S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
		Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
	}
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
	for(j in 1:(Nsamples*Nchains/2)){
		if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
			ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
			ll_mx[k] = max(abs(Lambda_new[,k,j]))
	}
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		for(k in 1:d){
			Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
		}
		Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
			Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
	}
}

## Varimax Rotation
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		tmp = varimax(Lambda_new[,,ind])
		Rot_vmx[,,ind] = tmp$rot
		Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
		Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
	}
}

## Permute factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
		Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
		Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
		}
	}
}

## Variance explained by retained factors
var_new = array(0, dim=c(d,Nsamples*Nchains/2))
var_vmx = array(0, dim=c(d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
		var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
	}
}

## Communalities & Uniquenesses
comm_uniq = array(0, dim=c(P,2,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		comm_uniq[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
	}
}

## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
	for(p2 in 1:P){
		Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
	}
}
comm_uniq_m = array(0,dim=c(P,2,3))
Lambda_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
	comm_uniq_m[p,1,] = quantile(comm_uniq[p,1,],probs=c(.025,.5,.975))
	comm_uniq_m[p,2,] = quantile(comm_uniq[p,2,],probs=c(.025,.5,.975))
	for(k in 1:d){
		Lambda_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
		Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
	}
}
Theta_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
	for(i in 1:N){
		Theta_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
		Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
	}
}		
var_splithalf2_new_m = array(0,dim=c(d,3))
var_splithalf2_vmx_m = array(0,dim=c(d,3))
for(k in 1:d){
	var_splithalf2_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
	var_splithalf2_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
}
var_splithalf2_expl = apply(var_splithalf2_new_m,2,sum)/P
var_splithalf2_vmx_expl = apply(var_splithalf2_vmx_m,2,sum)/P


lambda_splithalf2.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda_splithalf2.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda_splithalf2.tab) = c("",names(data)[ind_Y])
lambda_splithalf2.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda_splithalf2.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_vmx_m[,d1,1],3)
	lambda_splithalf2.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_vmx_m[,d1,2],3)
	lambda_splithalf2.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_vmx_m[,d1,3],3)
}

## Create SM Table 6
tabAsStringMatrix = print(lambda_splithalf2.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm6 = knitr::kable(tabAsStringMatrix)


########################################
########################################
## Impute Missing NP Data & Run BPPCA ## 
########################################
########################################

# Number of multiple imputed datasets & maximum number of iterations 
n.imp = 1
n.iter = 5

var.ls <- c("RaceEthnicity","HighestParentalEducation","HouseholdMaritalStatus","HouseholdIncome","PicVocab","Flanker",
				"List","CardSort","Pattern","Picture","Reading","RAVLT","LMT")
dat0 = data.table(data)
dat0 <- dat0[, var.ls, with = FALSE ]
Hmisc::describe(dat0)

# Initialize model
ini = mice( dat0, m = 1, maxit = 0 )
meth = ini$meth
meth["RaceEthnicity"] = "polyreg"
meth["HighestParentalEducation"] = "polyreg"
meth["HouseholdMaritalStatus"] = "polyreg"
meth["HouseholdIncome"] = "polyreg"
meth["PicVocab"] = "norm.predict"
meth["Flanker"] = "norm.predict"
meth["List"] = "norm.predict"
meth["CardSort"] = "norm.predict"
meth["Pattern"] = "norm.predict"
meth["Picture"] = "norm.predict"
meth["Reading"] = "norm.predict"
meth["RAVLT"] = "norm.predict"
meth["LMT"] = "norm.predict"
pred = ini$pred

set.seed(314)
post = mice( dat0, meth = meth, pred = pred, seed = 111,
              m = 1, maxit = 0)$post
dat.imp = mice( dat0, meth = meth, pred = pred, post = post,
                 seed = 1111, m = n.imp, maxit = n.iter)
rm(dat0)
# Convert imputed data to long format
dat.mi <- complete(dat.imp, action = "long", include = TRUE)
rm(dat.imp)

names(dat.mi)[1:2] <- c("imp", "id")
dat.mi$pid = rep(data$pid,n.imp+1)
dat.mi = dat.mi[dat.mi$imp==1,]
dat.mi = merge(dat.mi,data[,c(1,which(names(data)=="Female"),which(names(data)=="Age"),
			which(names(data)=="site_num"),which(names(data)=="rel_family_id"))],by="pid")
names(dat.mi)
dim(dat.mi)

###############
## Prep data ##
###############

ind_Y = c(8:16); names(dat.mi)[ind_Y]
ind_demog = c(4:7,17:18); names(dat.mi)[ind_demog]
ind_nest = c(19,20); names(dat.mi)[ind_nest]

data.mi1 = dat.mi[complete.cases(dat.mi[,c(ind_Y,ind_nest)]),]; names(data.mi1); dim(data.mi1)
site_num = rep(NA,length(data.mi1$site_num))
for(s in 1:length(unique(data.mi1$site_num))){
	site_num[data.mi1$site_num == unique(data.mi1$site_num)[s]] = s
}
data.mi1$site_num = site_num; rm(site_num)
data.mi1$id_fam = 0
data.mi1$fam_size = 0
ind=0
for(s in 1:length(unique(data.mi1$site_num))){
	data_s = data.mi1[data.mi1$site_num == s, ]
	for(f in 1:length(unique(data_s$rel_family_id))){
		data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
			sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
		if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
			ind=ind+1	
			data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
		}	
	}
	data.mi1[data.mi1$site_num == s, ] = data_s
}
data.mi1 = data.mi1[order(data.mi1$site_num,data.mi1$id_fam),]

Site = data.mi1$site_num
Fam = data.mi1$id_fam
ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
	ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
	nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)
Y = (as.matrix(scale(data.mi1[,8:16]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

####################################
## Run stan model on imputed data ##
####################################

Nsamples = 1000
Nchains = 2

set.seed(314)
model_file = "stan_code/bppca.stan"
smod = stan_model(model_file)

D_max = 3
sa.imputed.list = list()
log_lik.imputed.list = list()
looic.imputed.list = list()
for(d in D_max){
	print(d)
	pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
	set.seed(314)
	sa.imputed.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")  
	log_lik.imputed.list[[d]] <- extract(sa.imputed.list[[d]],"log_lik_marg")[[1]]
	looic.imputed.list[[d]] = loo(log_lik.imputed.list[[d]])
	save(sa.imputed.list,log_lik.imputed.list,
		looic.imputed.list,file="results/bppca_imputed_results.RData")
	print("###############################")
}

#looic.imputed.obj = compare(looic.imputed.list[[1]],looic.imputed.list[[2]],looic.imputed.list[[3]],
#	looic.imputed.list[[4]],looic.imputed.list[[5]],looic.imputed.list[[6]],
#	looic.imputed.list[[7]],looic.imputed.list[[8]])
#print(looic.imputed.obj)

###################################
## Extract parameters from model ##
###################################

load("results/bppca_imputed_results.RData")

d = 3
sa = sa.imputed.list[[d]]

Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

print(sa, pars=c("Q"), probs=c(.025,.5,.975))
print(sa, pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"), probs=c(.025,.5,.975))
traceplot(sa, pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"))

vc_imputed_tab = array(0, dim=c(4,3))
vc_imputed_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_imputed_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_imputed_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_imputed_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nsamples*Nchains/2))
R_new = array(0, dim=c(P,P,Nsamples*Nchains/2))
S_new = array(0, dim=c(d,d,Nsamples*Nchains/2)) 
W_new = array(0, dim=c(d,d,Nsamples*Nchains/2)) 
Theta_old = array(0, dim=c(d,N,Nsamples*Nchains/2))
Theta_new = array(0, dim=c(d,N,Nsamples*Nchains/2))
Lambda_old = array(0, dim=c(P,d,Nsamples*Nchains/2))
Lambda_new = array(0, dim=c(P,d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind+ 1
		Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
		Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
		Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
		W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
		R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
		S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
		Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
	}
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
	for(j in 1:(Nsamples*Nchains/2)){
		if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
			ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
			ll_mx[k] = max(abs(Lambda_new[,k,j]))
	}
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		for(k in 1:d){
			Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
		}
		Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
			Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
	}
}

## Varimax Rotation
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		tmp = varimax(Lambda_new[,,ind])
		Rot_vmx[,,ind] = tmp$rot
		Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
		Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
	}
}

## Permute factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
		Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
		Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
		}
	}
}

## Variance explained by retained factors
var_new = array(0, dim=c(d,Nsamples*Nchains/2))
var_vmx = array(0, dim=c(d,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
		var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
	}
}

## Communalities & Uniquenesses
comm_uniq = array(0, dim=c(P,2,Nsamples*Nchains/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		comm_uniq[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
	}
}

## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
	for(p2 in 1:P){
		Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
	}
}
comm_uniq_m = array(0,dim=c(P,2,3))
Lambda_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
	comm_uniq_m[p,1,] = quantile(comm_uniq[p,1,],probs=c(.025,.5,.975))
	comm_uniq_m[p,2,] = quantile(comm_uniq[p,2,],probs=c(.025,.5,.975))
	for(k in 1:d){
		Lambda_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
		Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
	}
}
Theta_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
	for(i in 1:N){
		Theta_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
		Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
	}
}		
var_imputed_new_m = array(0,dim=c(d,3))
var_imputed_vmx_m = array(0,dim=c(d,3))
for(k in 1:d){
	var_imputed_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
	var_imputed_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
}
var_imputed_expl = apply(var_imputed_new_m,2,sum)/P
var_imputed_vmx_expl = apply(var_imputed_vmx_m,2,sum)/P

lambda_imputed.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda_imputed.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda_imputed.tab) = c("",names(data.mi1)[ind_Y])
lambda_imputed.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda_imputed.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_vmx_m[,d1,1],3)
	lambda_imputed.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_vmx_m[,d1,2],3)
	lambda_imputed.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_vmx_m[,d1,3],3)
}

## Create SM Table 7
tabAsStringMatrix = print(lambda_imputed.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm7 = knitr::kable(tabAsStringMatrix)



##############################
##############################
## Export tables for paper  ##
##############################
##############################

save(file="figures_and_tables/bppca_tables_final.RData", 
	tab1,tab2,tab3,tab4,tab5,tab6,tab_sm1a,tab_sm1b,tab_sm2a,tab_sm2b,
	tab_sm3,tab_sm4,tab_sm5,tab_sm6,tab_sm7,looic.obj,vc_tab,
	var_expl_new,var_expl_vmx,var_expl_pmx,var_new_m,var_vmx_m,var_pmx_m,
	vc_nihtb_tab,var_nihtb_expl_new,var_nihtb_expl_vmx,var_nihtb_expl_pmx,
	var_nihtb_new_m,var_nihtb_vmx_m,var_nihtb_pmx_m,looic.nihtb.obj,
	d_rsq,d_rsq_nihtbx,d_rsq.cv,d_rsq.cv.nihtbx,beta.boot)

