rm(list = ls())


install.packages("Rtools")
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("forestplot")
install.packages("forestploter")
install.packages("MVMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
install.packages("openxlsx")
install.packages("ieugwasr")

devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()
library(openxlsx)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(forestplot)
library(forestploter)
library(metafor)
library(RadialMR)
library(MRPRESSO)
library(MVMR)
library(ieugwasr)

setwd("UKB/results")

###VAT combined########### 


###read in exposure data
vat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/VAT_GWAS_effects.csv')
exp_VAT <- format_data(vat_exp, type = "exposure")

##read in outcome data

vat_gs_out<-read.csv('input_files/out_vat_grip_strength_sex_combined.csv')
out_VAT_gs <- format_data(vat_gs_out, type = "outcome")

##harmonise data

dat <- harmonise_data(
  exposure_dat = exp_VAT, 
  outcome_dat = out_VAT_gs, action=1
)


##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"


##write harmonised SNPs and results to file
write.xlsx(dat, "output_files/VAT_MR_harmonised_sex_combined.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/VAT_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.008)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))


##leave one out analysis

loo_VAT_c<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_VAT_c)
p3[[1]]

write.csv(loo_VAT_c, file = "output_files/loo_VAT.csv", row.names = FALSE)


###VAT females########### 


###read in exposure data
vat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/VAT_GWAS_effects.csv')
exp_VAT <- format_data(vat_exp, type = "exposure")

##read in outcome data

vat_gs_out<-read.csv('input_files/out_vat_grip_strength.csv')
vat_gs_out_female<-subset(vat_gs_out,sex==0)

out_VAT_gs_female <- format_data(vat_gs_out_female, type = "outcome")

##harmonise data

dat <- harmonise_data(
  exposure_dat = exp_VAT, 
  outcome_dat = out_VAT_gs_female, action=1
)


##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"

##Write results files 
write.xlsx(dat, "output_files/VAT_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/VAT_MR_results_female.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)


###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.008)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

loo_VAT_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_VAT_f)
p3[[1]]
write.csv(loo_VAT_f, file = "output_files/loo_VAT_f.csv", row.names = FALSE)

#### VAT males####

## read in exposure data
vat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/VAT_GWAS_effects.csv')
exp_VAT <- format_data(vat_exp, type = "exposure")

##read in outcome data

vat_gs_out<-read.csv('input_files/out_vat_grip_strength.csv')
vat_gs_out_male<-subset(vat_gs_out,sex==1)
out_VAT_gs_male <- format_data(vat_gs_out_male, type = "outcome")

##harmonise data

dat <- harmonise_data(
  exposure_dat = exp_VAT, 
  outcome_dat = out_VAT_gs_male, action=1
)

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##write results files
write.xlsx(dat, "output_files/VAT_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/VAT_MR_results_male.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)


###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.008)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

loo_VAT_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_VAT_m)
p3[[1]]
write.csv(loo_VAT_m, file = "output_files/loo_VAT_m.csv", row.names = FALSE)



###ASAT sex combined ######

###read in exposure data

asat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/ASAT_GWAS_effects.csv')
exp_ASAT <- format_data(asat_exp, type = "exposure")

##read in outcome data

asat_gs_out<-read.csv('input_files/out_asat_grip_strength_sex_combined.csv')
out_ASAT_gs <- format_data(asat_gs_out, type = "outcome")

##harmonise data

dat <- harmonise_data(
  exposure_dat = exp_ASAT, 
  outcome_dat = out_ASAT_gs, action=1
)

##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##write results files 
write.xlsx(dat, "output_files/ASAT_MR_harmonised_sex_combined.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/ASAT_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


##leave one out analysis
loo_ASAT_c<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_ASAT_c)
p3[[1]]
write.csv(loo_ASAT_c, file = "output_files/loo_ASAT.csv", row.names = FALSE)

###ASAT females ######


###read in exposure data

asat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/ASAT_GWAS_effects.csv')
exp_ASAT <- format_data(asat_exp, type = "exposure")

##read in outcome data

asat_gs_out<-read.csv('input_files/out_asat_grip_strength.csv')
asat_gs_out_female<-subset(asat_gs_out,sex==0)
out_ASAT_gs_female <- format_data(asat_gs_out_female, type = "outcome")


##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ASAT, 
  outcome_dat = out_ASAT_gs_female, action=1
)

## perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"

##write results files
write.xlsx(dat, "output_files/ASAT_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/ASAT_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


##leave one out analysis
loo_ASAT_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_ASAT_f)
p3[[1]]
write.csv(loo_ASAT_f, file = "output_files/loo_ASAT_f.csv", row.names = FALSE)

#####ASAT males############################

###read in exposure data

asat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/ASAT_GWAS_effects.csv')
exp_ASAT <- format_data(asat_exp, type = "exposure")


##read in outcome data
asat_gs_out<-read.csv('input_files/out_asat_grip_strength.csv')
asat_gs_out_male<-subset(asat_gs_out,sex==1)
out_ASAT_gs_male <- format_data(asat_gs_out_male, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ASAT, 
  outcome_dat = out_ASAT_gs_male, action=1
)

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##write results files
write.xlsx(dat, "output_files/ASAT_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/ASAT_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


##leave one out analysis
loo_ASAT_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_ASAT_m)
p3[[1]]
write.csv(loo_ASAT_m, file = "output_files/loo_ASAT_m.csv", row.names = FALSE)

#####GFAT combined#####

###read in exposure data
gfat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/GFAT_GWAS_effects.csv')
exp_GFAT<- format_data(gfat_exp, type = "exposure")


##read in outcome data
gfat_gs_out<-read.csv('input_files/out_gfat_grip_strength_sex_combined.csv')
out_GFAT_gs<- format_data(gfat_gs_out, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_GFAT, 
  outcome_dat = out_GFAT_gs, action=1
)

##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_ivw_mre", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##write results files
write.xlsx(dat, "output_files/GFAT_MR_harmonised.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/GFAT_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)



##leave one out analysis

loo_GFAT_c<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_GFAT_c)
p3[[1]]
write.csv(loo_GFAT_c, file = "output_files/loo_GFAT_c.csv", row.names = FALSE)


#####GFAT females#####

###read in exposure data
gfat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/GFAT_GWAS_effects.csv')
exp_GFAT<- format_data(gfat_exp, type = "exposure")


##read in outcome data

gfat_gs_out<-read.csv('input_files/out_gfat_grip_strength.csv')
gfat_gs_out_female<-subset(gfat_gs_out,sex==0)
out_GFAT_gs_female <- format_data(gfat_gs_out_female, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_GFAT, 
  outcome_dat = out_GFAT_gs_female, action=1
)

##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"


##write results to file
write.xlsx(dat, "output_files/GFAT_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/GFAT_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

loo_GFAT_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_GFAT_f)
p3[[1]]
write.csv(loo_GFAT_f, file = "output_files/loo_GFAT_females.csv", row.names = FALSE)



#####GFAT males ####

###read in exposure data
gfat_exp<-read.csv('GWAS/VAT_ASAT_GFAT/GFAT_GWAS_effects.csv')
exp_GFAT<- format_data(gfat_exp, type = "exposure")


##read in outcome data
gfat_gs_out<-read.csv('input_files/out_gfat_grip_strength.csv')
gfat_gs_out_male<-subset(gfat_gs_out,sex==1)
out_GFAT_gs_male <- format_data(gfat_gs_out_male, type = "outcome")


##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_GFAT, 
  outcome_dat = out_GFAT_gs_male, action=1
)
dat<-subset(dat, dat$SNP!="rs200472737")

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##
write.xlsx(dat, "output_files/GFAT_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/GFAT_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

loo_GFAT_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_GFAT_m)
p3[[1]]
write.csv(loo_GFAT_m, file = "output_files/loo_GFAT_males.csv", row.names = FALSE)



##BF combined#######

##read in exposure data
ukbf_exp<-read.csv('GWAS/Bodyfat/BF_GWAS_effects.csv')
exp_ukbf <- format_data(ukbf_exp, type = "exposure")


##read in outcome data
ukbf_gs_out<-read.csv('input_files/exp_ukbf_SNP_gs_sex_combined.csv')
out_ukbf_gs <- format_data(ukbf_gs_out, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ukbf, 
  outcome_dat = out_ukbf_gs, action=1
)

##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_ivw_mre", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##write results to file
write.xlsx(dat, "output_files/BF_MR_harmonised_sex_combined.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/BF_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.0002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

loo_BF_c<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_BF_c)
p3[[1]]
write.csv(loo_BF_c, file = "output_files/loo_BF_c.csv", row.names = FALSE)


##BF females#######

##read in exposure data
ukbf_exp<-read.csv('GWAS/Bodyfat/BF_GWAS_effects.csv')
exp_ukbf <- format_data(ukbf_exp, type = "exposure")


##read in outcome data
ukbf_gs_out<-read.csv('input_files/exp_ukbf_SNP_gs.csv')
ukbf_gs_out_female<-subset(ukbf_gs_out,sex==0)

out_ukbf_gs_female <- format_data(ukbf_gs_out_female, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ukbf, 
  outcome_dat = out_ukbf_gs_female, action=1
)

##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"


##write results files
write.xlsx(dat, "output_files/BF_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/BF_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.0002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)




##leave one out analysis

loo_BF_f<-mr_leaveoneout(dat)
write.csv(loo_BF_f, file = "output_files/loo_BF_females.csv", row.names = FALSE)

##BF males#######

ukbf_exp<-read.csv('C:/Users/rmjlaet/OneDrive - University College London/Snehal/Adiposity-strength-crf/GWAS/Bodyfat/BF_GWAS_effects.csv')
exp_ukbf <- format_data(ukbf_exp, type = "exposure")


##read in outcome data
ukbf_gs_out<-read.csv('input_files/exp_ukbf_SNP_gs.csv')
ukbf_gs_out_male<-subset(ukbf_gs_out,sex==1)
out_ukbf_gs_male <- format_data(ukbf_gs_out_male, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ukbf, 
  outcome_dat = out_ukbf_gs_male, action=1
)

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##write results files 
write.xlsx(dat, "output_files/BF_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/BF_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.0002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis
loo_BF_m<-mr_leaveoneout(dat)
write.csv(loo_BF_m, file = "output_files/loo_BF_males.csv", row.names = FALSE)


##FA sex combined#######

fa_exp<-read.csv('GWAS/Bodyfat/FA_GWAS_effects.csv')
exp_fa <- format_data(fa_exp, type = "exposure")


##read in outcome data
fa_gs_out<-read.csv('input_files/out_fa_grip_strength_sex_combined.csv')
out_fa_gs_comb <- format_data(fa_gs_out, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_fa, 
  outcome_dat = out_fa_gs_comb, action=1
)

##remove SNP in LD
dat<-subset(dat, dat$SNP!="rs12369179")

##perform MR analyses
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##write results to file
write.xlsx(dat, "output_files/FA_MR_harmonised.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/FA_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.001)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

loo_FA<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_FA)
p3[[1]]
write.csv(loo_FA, file = "output_files/loo_FA.csv", row.names = FALSE)


##FA females#######

##read exposure data
fa_exp<-read.csv('GWAS/Bodyfat/FA_GWAS_effects.csv')
exp_fa <- format_data(fa_exp, type = "exposure")


##read in outcome data
fa_gs_out<-read.csv('input_files/out_fa_grip_strength.csv')
fa_gs_out_female<-subset(fa_gs_out,sex==0)
out_fa_gs_female <- format_data(fa_gs_out_female, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_fa, 
  outcome_dat = out_fa_gs_female, action=1
)

##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"

##write results files
write.xlsx(dat, "output_files/FA_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/FA_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.001)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

loo_FA_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_FA_f)
p3[[1]]
write.csv(loo_FA_f, file = "output_files/loo_FA_females.csv", row.names = FALSE)


##FA males#######

##read exposure data
fa_exp<-read.csv('GWAS/Bodyfat/FA_GWAS_effects.csv')
exp_fa <- format_data(fa_exp, type = "exposure")


##read in outcome data

fa_gs_out<-read.csv('input_files/out_fa_grip_strength.csv')
#fa_gs_out_female<-subset(fa_gs_out,sex==0)
fa_gs_out_male<-subset(fa_gs_out,sex==1)

#out_fa_gs_female <- format_data(fa_gs_out_female, type = "outcome")
out_fa_gs_male <- format_data(fa_gs_out_male, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_fa, 
  outcome_dat = out_fa_gs_male, action=1
)

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##write results files
write.xlsx(dat, "output_files/FA_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/FA_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.001)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

loo_FA_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_FA_m)
p3[[1]]
write.csv(loo_FA_m, file = "output_files/loo_FA_males.csv", row.names = FALSE)





##UFA sex combined#######

##read in exposure data
ufa_exp<-read.csv('GWAS/Bodyfat/UFA_GWAS_effects.csv')
exp_ufa <- format_data(ufa_exp, type = "exposure")


##read in outcome data
ufa_gs_out<-read.csv('input_files/out_ufa_grip_strength_sex_combined.csv')
out_ufa_gs_comb <- format_data(ufa_gs_out, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ufa, 
  outcome_dat = out_ufa_gs_comb, action=1
)

##subset data to remove SNP in LD
dat <-subset(dat, dat$SNP!="chr1_72767554")

##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##Write results to file
write.xlsx(dat, "output_files/UFA_MR_harmonised.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/UFA_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.001)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

dat$exposure<-"MetUFA"

loo_UFA<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_UFA)
p3[[1]]
write.csv(loo_UFA, file = "output_files/loo_UFA_c.csv", row.names = FALSE)

##UFA females#######

ufa_exp<-read.csv('GWAS/Bodyfat/UFA_GWAS_effects.csv')
exp_ufa <- format_data(ufa_exp, type = "exposure")

##read in outcome data
ufa_gs_out<-read.csv('input_files/out_ufa_grip_strength.csv')
ufa_gs_out_female<-subset(ufa_gs_out,sex==0)
out_ufa_gs_female <- format_data(ufa_gs_out_female, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ufa, 
  outcome_dat = out_ufa_gs_female, action=1
)

##subset to remove SNP in LD
dat <-subset(dat, dat$SNP!="chr1_72767554")

##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"

##write results to file
write.xlsx(dat, "output_files/UFA_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/UFA_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.001)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis
dat$exposure<-"MetUFA"

loo_UFA_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_UFA_f)
p3[[1]]
write.csv(loo_UFA_f, file = "output_files/loo_UFA_females.csv", row.names = FALSE)


##UFA males#######

##read in exposure data
ufa_exp<-read.csv('GWAS/Bodyfat/UFA_GWAS_effects.csv')
exp_ufa <- format_data(ufa_exp, type = "exposure")


##read in outcome data
ufa_gs_out<-read.csv('input_files/out_ufa_grip_strength.csv')
ufa_gs_out_male<-subset(ufa_gs_out,sex==1)
out_ufa_gs_male <- format_data(ufa_gs_out_male, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ufa, 
  outcome_dat = out_ufa_gs_male, action=1
)

##remove SNP in LD
dat <-subset(dat, dat$SNP!="chr1_72767554")

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##write results to file
write.xlsx(dat, "output_files/UFA_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/UFA_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.001)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

dat$exposure<-"MetUFA"

loo_UFA_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_UFA_m)
p3[[1]]
write.csv(loo_UFA_m, file = "output_files/loo_UFA_males.csv", row.names = FALSE)

#####ATMFI combined#####

###read in exposure data
ATMFI_exp<-read.csv('GWAS/ATFMI/ATMFI_GWAS_effects.csv')
ATMFI_exp$Phenotype<-"ATMFI"
exp_ATMFI<- format_data(ATMFI_exp, type = "exposure")

##read in outcome data
atmfi_gs_out<-read.csv('input_files/out_ATMFI_grip_strength_sex_combined.csv')
out_ATMFI_gs <- format_data(atmfi_gs_out, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ATMFI, 
  outcome_dat = out_ATMFI_gs, action=1
)



##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##write results files
write.xlsx(dat, "output_files/ATMFI_MR_harmonised_sex_combined.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/ATMFI_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.003)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)

##leave one out analysis

loo_ATMFI_c<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_ATMFI_c)
p3[[1]]
write.csv(loo_ATMFI_c, file = "output_files/loo_ATMFI.csv", row.names = FALSE)

#####ATMFI females#####

###read in exposure data
ATMFI_exp<-read.csv('GWAS/ATFMI/ATMFI_GWAS_effects.csv')
ATMFI_exp$Phenotype<-"ATMFI"
exp_ATMFI<- format_data(ATMFI_exp, type = "exposure")

##read in outcome data
atmfi_gs_out<-read.csv('input_files/out_ATMFI_grip_strength.csv')
atmfi_gs_out_female<-subset(atmfi_gs_out,sex==0)
out_ATMFI_gs_female <- format_data(atmfi_gs_out_female, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ATMFI, 
  outcome_dat = out_ATMFI_gs_female, action=1
)

##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"

##write results to file
write.xlsx(dat, "output_files/ATMFI_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/ATMFI_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.003)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis
loo_ATMFI_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_ATMFI_f)
p3[[1]]
write.csv(loo_ATMFI_f, file = "output_files/loo_ATMFI_females.csv", row.names = FALSE)

####ATMFI males#####

###read in exposure data
ATMFI_exp<-read.csv('GWAS/ATFMI/ATMFI_GWAS_effects.csv')
ATMFI_exp$Phenotype<-"ATMFI"
exp_ATMFI<- format_data(ATMFI_exp, type = "exposure")


##read in outcome data
atmfi_gs_out<-read.csv('input_files/out_ATMFI_grip_strength.csv')
atmfi_gs_out_male<-subset(atmfi_gs_out,sex==1)
out_ATMFI_gs_male <- format_data(atmfi_gs_out_male, type = "outcome")

ATMFI_exp$Phenotype<-"ATMFI"

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_ATMFI, 
  outcome_dat = out_ATMFI_gs_male, action=1
)

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

##write results to file
write.xlsx(dat, "output_files/ATMFI_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"M"
write.table(res_males, "output_files/ATMFI_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.003)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis
loo_ATMFI_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_ATMFI_m)
p3[[1]]
write.csv(loo_ATMFI_m, file = "output_files/loo_ATMFI_males.csv", row.names = FALSE)


#####PTMFI sex combined#####

###read in exposure data
ptmfi_exp<-read.csv('GWAS/PTFMI/PTMFI_GWAS_effects.csv')
exp_PTMFI<- format_data(ptmfi_exp, type = "exposure")


##read in outcome data
ptmfi_gs_out<-read.csv('input_files/out_PTMFI_grip_strength_sex_combined.csv')
out_PTMFI_gs <- format_data(ptmfi_gs_out, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_PTMFI, 
  outcome_dat = out_PTMFI_gs, action=1
)



dat$exposure<-"PTMFI"
dat$outcome<-"Grip Strength"

##perform MR analysis
res_comb <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"Both"

##write results to file
write.xlsx(dat, "output_files/PTMFI_MR_harmonised_sex_combined.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_comb$sex<-"Both"
write.table(res_comb, "output_files/PTMFI_MR_results.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

loo_PTMFI_c<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_PTMFI_c)
p3[[1]]
write.csv(loo_PTMFI_c, file = "output_files/loo_PTMFI.csv", row.names = FALSE)

#####PTMFI females#####

###read in exposure data
ptmfi_exp<-read.csv('GWAS/PTFMI/PTMFI_GWAS_effects.csv')
exp_PTMFI<- format_data(ptmfi_exp, type = "exposure")


##read in outcome data

ptmfi_gs_out<-read.csv('input_files/out_PTMFI_grip_strength.csv')
ptmfi_gs_out_female<-subset(ptmfi_gs_out,sex==0)

out_PTMFI_gs_female <- format_data(ptmfi_gs_out_female, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_PTMFI, 
  outcome_dat = out_PTMFI_gs_female, action=1
)

##perform MR analysis
res_females <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"F"

dat$exposure<-"PTMFI"
dat$outcome<-"Grip Strength"

##write results to file
write.xlsx(dat, "output_files/PTMFI_MR_harmonised_females.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_females$sex<-"F"
write.table(res_females, "output_files/PTMFI_MR_results_females.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

loo_PTMFI_f<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_PTMFI_f)
p3[[1]]
write.csv(loo_PTMFI_f, file = "output_files/loo_PTMFI_females.csv", row.names = FALSE)

#####PTMFI males#####

###read in exposure data
ptmfi_exp<-read.csv('GWAS/PTFMI/PTMFI_GWAS_effects.csv')
exp_PTMFI<- format_data(ptmfi_exp, type = "exposure")


##read in outcome data
ptmfi_gs_out<-read.csv('input_files/out_PTMFI_grip_strength.csv')
ptmfi_gs_out_male<-subset(ptmfi_gs_out,sex==1)
out_PTMFI_gs_male <- format_data(ptmfi_gs_out_male, type = "outcome")

##harmonise data
dat <- harmonise_data(
  exposure_dat = exp_PTMFI, 
  outcome_dat = out_PTMFI_gs_male, action=1
)

##perform MR analysis
res_males <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

##calculate F-statistics
dat$f<-(dat$beta.exposure/dat$se.exposure)^2
dat$meanf<-mean(dat$f)
dat$sex<-"M"

dat$exposure<-"PTMFI"
dat$outcome<-"Grip Strength"

##write results to file
write.xlsx(dat, "output_files/PTMFI_MR_harmonised_males.xlsx", sep = ",", colNames = TRUE, rowNames = FALSE, append = TRUE)
res_males$sex<-"M"
write.table(res_males, "output_files/PTMFI_MR_results_males.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

###heterogeneity
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

###radial MR 
radial_dat<-format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
radial_res<-ivw_radial(radial_dat, 0.002)
snp_list_radial_res <- radial_res[["outliers"]][["SNP"]]
dat_filtered <- dat[!dat$SNP %in% snp_list_radial_res, ]

res_rad_exc<- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
mr_heterogeneity(dat_filtered)


##leave one out analysis

loo_PTMFI_m<-mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(loo_PTMFI_m)
p3[[1]]
write.csv(loo_PTMFI_m, file = "output_files/loo_PTMFI_males.csv", row.names = FALSE)