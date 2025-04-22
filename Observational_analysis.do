
**SD of adiposity measures- sex combined- adjusted 

cd "$od"

tempname vat

postfile `vat' str15 Adiposity N beta se using Observational_SD_sex_combined.dta, replace
foreach var of varlist VAT_sd ASAT_sd BF_sd2 GFAT_sd ATMFI_sd PTMFI_sd {
	regress gs_sd2 `var' i.smoking_status2 i.townsend_5 i.ls_pa_4days_iv i.ls_alc_freq3cats_iv ls_vpa_4days_iv age2 sex height if eur_unrelated==1 & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var', N: `N'"
	post  `vat'  ("`var'") (`N') (_b[`var']) (_se[`var'])
}
 

postclose `vat'

**SD of adiposity measures- sex combined- unadjusted 

cd "$od"

tempname vat

postfile `vat' str15 Adiposity N beta se using Observational_SD_sex_combined_unadjusted.dta, replace
foreach var of varlist VAT_sd ASAT_sd BF_sd2 GFAT_sd ATMFI_sd PTMFI_sd {
	regress gs_sd2 `var' if eur_unrelated==1 & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var', N: `N'"
	post  `vat'  ("`var'") (`N') (_b[`var']) (_se[`var'])
}
 

postclose `vat'



**SD of adiposity measures- sex stratified- adjusted

cd "$od"

tempname vat

postfile `vat' str15 Adiposity sex N beta se using Observational_SD.dta, replace
foreach var of varlist VAT_sd ASAT_sd BF_sd2 GFAT_sd ATMFI_sd PTMFI_sd {
	forvalues sex=0/1{
	regress gs_sd2 `var' i.smoking_status2 i.townsend_5 i.ls_pa_4days_iv i.ls_alc_freq3cats_iv ls_vpa_4days_iv age2 height if eur_unrelated==1 & sex==`sex' & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var', Sex: `sex', N: `N'"
	post  `vat'  ("`var'") (`sex') (`N') (_b[`var']) (_se[`var'])
}
 
}

postclose `vat'

**SD of adiposity measures- sex stratified- unadjusted

cd "$od"

tempname vat

postfile `vat' str15 Adiposity sex N beta se using Observational_SD_unadj.dta, replace
foreach var of varlist VAT_sd ASAT_sd BF_sd2 GFAT_sd ATMFI_sd PTMFI_sd {
	forvalues sex=0/1{
	regress gs_sd2 `var' if eur_unrelated==1 & sex==`sex' & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var', Sex: `sex', N: `N'"
	post  `vat'  ("`var'") (`sex') (`N') (_b[`var']) (_se[`var'])
}
 
}

postclose `vat'



**Quintiles of adiposity measures- sex combined- adjusted

cd "$od"

tempname vat

postfile `vat' str15 Adiposity N beta se using Observational_Quintiles_sex_combined.dta, replace
foreach var of varlist VAT_5 ASAT_5 bf_52 GFAT_5 ATMFI_5 PTMFI_5  {
	regress gs_sd2 i.`var' i.smoking_status2 i.townsend_5 i.ls_pa_4days_iv ls_vpa_4days_iv i.ls_alc_freq3cats_iv age2 sex height if eur_unrelated==1 & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var',  N: `N'"
	post  `vat'  ("`var'")  (`N') (_b[2.`var']) (_se[2.`var']) 
	post  `vat'  ("`var'")  (`N') (_b[3.`var']) (_se[3.`var'])
	post  `vat'  ("`var'") (`N') (_b[4.`var']) (_se[4.`var']) 
	post  `vat'  ("`var'")  (`N') (_b[5.`var']) (_se[5.`var']) 
}
 

postclose `vat'

**Quintiles of adiposity measures- sex combined- unadjusted

cd "$od"

tempname vat

postfile `vat' str15 Adiposity N beta se using Observational_Quintiles_sex_combined_unadjusted.dta, replace
foreach var of varlist VAT_5 ASAT_5 bf_52 GFAT_5 ATMFI_5 PTMFI_5  {
	regress gs_sd2 i.`var' if eur_unrelated==1 & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var',  N: `N'"
	post  `vat'  ("`var'")  (`N') (_b[2.`var']) (_se[2.`var']) 
	post  `vat'  ("`var'")  (`N') (_b[3.`var']) (_se[3.`var'])
	post  `vat'  ("`var'") (`N') (_b[4.`var']) (_se[4.`var']) 
	post  `vat'  ("`var'")  (`N') (_b[5.`var']) (_se[5.`var']) 
}
 

postclose `vat'


**Quintiles of adiposity measures- sex stratified- adjusted

cd "$od"

tempname vat

postfile `vat' str15 Adiposity sex N beta se using Observational_Quintiles.dta, replace
foreach var of varlist VAT_5 ASAT_5 bf_52 GFAT_5 ATMFI_5 PTMFI_5 {
	forvalues sex=0/1{
	regress gs_sd2 i.`var' i.smoking_status2 i.townsend_5 i.ls_pa_4days_iv ls_vpa_4days_iv i.ls_alc_freq3cats_iv age2 height if eur_unrelated==1 & sex==`sex' & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var', Sex: `sex', N: `N'"
	post  `vat'  ("`var'")  (`sex') (`N') (_b[2.`var']) (_se[2.`var']) 
	post  `vat'  ("`var'") (`sex')  (`N') (_b[3.`var']) (_se[3.`var'])
	post  `vat'  ("`var'") (`sex') (`N') (_b[4.`var']) (_se[4.`var']) 
	post  `vat'  ("`var'") (`sex')  (`N') (_b[5.`var']) (_se[5.`var']) 
}
 
}

postclose `vat'

**Quintiles of adiposity measures- sex stratified- unadjusted

cd "$od"

tempname vat

postfile `vat' str15 Adiposity sex N beta se using Observational_Quintiles_unadj.dta, replace
foreach var of varlist VAT_5 ASAT_5 bf_52 GFAT_5 ATMFI_5 PTMFI_5  {
	forvalues sex=0/1{
	regress gs_sd2 i.`var' if eur_unrelated==1 & sex==`sex' & any_MRI_data==1
	local N=e(N) 
	display "Variable: `var', Sex: `sex', N: `N'"
	post  `vat'  ("`var'")  (`sex') (`N') (_b[2.`var']) (_se[2.`var']) 
	post  `vat'  ("`var'") (`sex')  (`N') (_b[3.`var']) (_se[3.`var'])
	post  `vat'  ("`var'") (`sex') (`N') (_b[4.`var']) (_se[4.`var']) 
	post  `vat'  ("`var'") (`sex')  (`N') (_b[5.`var']) (_se[5.`var']) 
}
 
}

postclose `vat'





