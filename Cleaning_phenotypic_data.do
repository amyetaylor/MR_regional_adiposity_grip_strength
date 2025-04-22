***This Stata do file cleans the UKB phenotypic data 

***sex*******************************************

label define sexlab 0 "Female" 1"Male"
label values n_31_0_0 sexlab
rename n_31_0_0 sex


****age*****************************************
**age at baseline assessment
rename n_21022_0_0 age

**age at imaging assessment visit
gen age2= n_21003_2_0

**age categories
gen agecat2= age2
recode agecat2 min/55=1 56/65=2 66/75=3 76/max=4
label define agecatlab 1"55 and under" 2"56-65" 3"66-75" 4"76 and over"
label values agecat agecatlab

***townsend (only one measure at baseline)***************

gen townsend=n_22189_0_0

xtile townsend_5=n_22189_0_0, nq(5)
lab def sep_townsendquint_bllabel 1"Least deprived" 2"2nd least deprived" 3"Median deprivation level" 4"2nd most deprived" 5"Most deprived"
label values townsend_5 sep_townsendquint_bllabel

***height 

gen height= n_50_2_0/100
replace height= n_50_0_0/100 if height==.


****grip strength**********************************************************************

**grip strength at baseline
gen grip_strength= n_46_0_0 if n_46_0_0!=. & n_46_0_0!=0
replace grip_strength = n_47_0_0 if n_46_0_0<=n_47_0_0 & n_47_0_0!=.  
replace grip_strength= n_47_0_0 if grip_strength==. 
replace grip_strength=. if grip_strength==0

egen gs_sdf=std(grip_strength) if sex==0
egen gs_sdm=std(grip_strength) if sex==1
gen gs_sd= gs_sdf
replace gs_sd=gs_sdm if sex==1

**grip strength at imaging visit
gen grip_strength2= n_46_2_0 if n_46_2_0!=. & n_46_2_0!=0
replace grip_strength2 = n_47_2_0 if n_46_2_0<=n_47_2_0 & n_47_2_0!=.  
replace grip_strength2= n_47_2_0 if grip_strength==. 
replace grip_strength2=. if grip_strength2==0

egen gs_sdf2=std(grip_strength2) if sex==0
egen gs_sdm2=std(grip_strength2) if sex==1
gen gs_sd2= gs_sdf2
replace gs_sd2=gs_sdm2 if sex==1

******BMI*******************************************************************************

**bmi categories
rename n_21001_0_0 bmi
gen bmicat=bmi
recode bmicat min/18.49999=1 18.500000/24.99999=2 25.000000/29.999999=3 30.00000/max=4
label define bmicat 1"Underweight" 2"Normal" 3"Overweight" 4"Obese"
label values bmicat bmicat

**3 levels of bmicat
gen bmicat3=bmicat
recode bmicat 1=2

**quintiles of BMI
xtile bmi_5f=bmi if sex==0, nq(5)
xtile bmi_5m=bmi if sex==1, nq(5)
gen bmi_5= bmi_5f 
replace bmi_5= bmi_5m if bmi_5==.

**bmi at imaging visit
gen bmi2= n_21001_2_0


***Adiposity measures
*****BF, VAT and ASAT***********************************************************************


***Bodyfat from bioimpedance
**bf is bf at baseline, bf2 is bf at imaging visit
gen bf= n_23099_0_0
gen bf2= n_23099_2_0

xtile bf_5f=bf if sex==0, nq(5)
xtile bf_5m=bf if sex==1, nq(5)
gen bf_5= bf_5f 
replace bf_5= bf_5m if bf_5==.

egen BF_sdf=std(bf) if sex==0
egen BF_sdm=std(bf) if sex==1
gen BF_sd= BF_sdf
replace BF_sd=BF_sdm if sex==1

xtile bf_5f2=bf2 if sex==0, nq(5)
xtile bf_5m2=bf2 if sex==1, nq(5)
gen bf_52= bf_5f2 
replace bf_52= bf_5m2 if bf_52==.

egen BF_sdf2=std(bf2) if sex==0
egen BF_sdm2=std(bf2) if sex==1
gen BF_sd2= BF_sdf2
replace BF_sd2=BF_sdm2 if sex==1

***VAT*******************************

gen VAT= n_22407_2_0
replace VAT=. if VAT_error!=.

xtile VAT_5f=VAT if sex==0, nq(5)
xtile VAT_5m=VAT if sex==1, nq(5)
gen VAT_5= VAT_5f 
replace VAT_5= VAT_5m if VAT_5==.

egen VAT_sdf=std(VAT) if sex==0 
egen VAT_sdm=std(VAT) if sex==1
gen VAT_sd= VAT_sdf
replace VAT_sd=VAT_sdm if sex==1

***ASAT******************************************

gen ASAT= n_22408_2_0
replace ASAT=. if ASAT_error!=.

egen ASAT_sdf=std(ASAT) if sex==0
egen ASAT_sdm=std(ASAT) if sex==1
gen ASAT_sd= ASAT_sdf
replace ASAT_sd=ASAT_sdm if sex==1

xtile ASAT_5f=ASAT if sex==0, nq(5)
xtile ASAT_5m=ASAT if sex==1, nq(5)
gen ASAT_5= ASAT_5f 
replace ASAT_5= ASAT_5m if ASAT_5==.

***ATMFI*************************************

gen ATMFI_l= AT_MFI_l if AT_error_l=="NA"
gen ATMFI_r= AT_MFI_r if AT_error_r=="NA"

gen ATMFI= (ATMFI_l + ATMFI_r)/2
replace ATMFI=ATMFI_l if ATMFI_r==.
replace ATMFI=ATMFI_r if ATMFI_l==.

egen ATMFI_sdf=std(ATMFI) if sex==0
egen ATMFI_sdm=std(ATMFI) if sex==1
gen ATMFI_sd= ATMFI_sdf
replace ATMFI_sd=ATMFI_sdm if sex==1

xtile ATMFI_5f=ATMFI if sex==0, nq(5)
xtile ATMFI_5m=ATMFI if sex==1, nq(5)
gen ATMFI_5= ATMFI_5f 
replace ATMFI_5= ATMFI_5m if ATMFI_5==.


**PTMFI****************************************

gen PTMFI_l= PT_MFI_l if PT_error_l=="NA"
gen PTMFI_r= PT_MFI_r if PT_error_r=="NA"

gen PTMFI= (PTMFI_l + PTMFI_r)/2
replace PTMFI=PTMFI_l if PTMFI_r==.
replace PTMFI=PTMFI_r if PTMFI_l==.

egen PTMFI_sdf=std(PTMFI) if sex==0
egen PTMFI_sdm=std(PTMFI) if sex==1
gen PTMFI_sd= PTMFI_sdf
replace PTMFI_sd=PTMFI_sdm if sex==1

xtile PTMFI_5f=PTMFI if sex==0, nq(5)
xtile PTMFI_5m=PTMFI if sex==1, nq(5)
gen PTMFI_5= PTMFI_5f 
replace PTMFI_5= PTMFI_5m if PTMFI_5==.

***GFAT***************************************

gen GFAT= TAT_volume- (VAT+ASAT)

egen GFAT_sdf=std(GFAT) if sex==0
egen GFAT_sdm=std(GFAT) if sex==1
gen GFAT_sd= GFAT_sdf
replace GFAT_sd=GFAT_sdm if sex==1

xtile GFAT_5f=GFAT if sex==0, nq(5)
xtile GFAT_5m=GFAT if sex==1, nq(5)
gen GFAT_5= GFAT_5f 
replace GFAT_5= GFAT_5m if GFAT_5==.



***Physical activity************************************************************************************

**moderate PA (days per week)
gen pa_mod=n_884_0_0 if  n_884_0_0>-1

**vigorous PA (days per week)
gen pa_vig=n_904_0_0 if  n_904_0_0>-1

*Physical activity 
 
**moderate physical activity- baseline
gen ls_pa_days_bl=n_884_0_0 
recode ls_pa_days_bl -7=. -1=. -3=.
lab var ls_pa_days_bl "Baseline number of days/week moderate physical activity >10 mins"
tab ls_pa_days_bl, mis
gen ls_pa_4days_bl=ls_pa_days_bl
recode ls_pa_4days_bl 7=0 6=0 5=0 4=0 3=1 2=1 1=1 0=1
lab var ls_pa_4days_bl "Baseline low PA - under 4 days/week moderate physical activity >10 mins"
lab def ls_pa_4days_bllab 0"Active - over 4 days/wk PA>10 mins" 1"Inactive - under 4 days/wk PA>10 mins"
lab val ls_pa_4days_bl ls_pa_4days_bllab

**vigorous physical activity- baseline
gen ls_vpa_days_bl=n_904_0_0 
recode ls_vpa_days_bl -7=. -1=. -3=.
lab var ls_vpa_days_bl "Baseline number of days/week vigorous physical activity >10 mins"
tab ls_vpa_days_bl, mis
gen ls_vpa_4days_bl=ls_vpa_days_bl
recode ls_vpa_4days_bl 7=0 6=0 5=0 4=0 3=1 2=1 1=1 0=1
lab var ls_vpa_4days_bl "Baseline low PA - under 4 days/week vigorous physical activity >10 mins"
lab def ls_vpa_4days_bllab 0"Active - over 4 days/wk PA>10 mins" 1"Inactive - under 4 days/wk PA>10 mins"
lab val ls_vpa_4days_bl ls_vpa_4days_bllab


**moderate physical activity- imaging visit
gen ls_pa_days_iv=n_884_2_0 
recode ls_pa_days_iv -7=. -1=. -3=.
lab var ls_pa_days_iv "Imaging visit number of days/week moderate physical activity >10 mins"
tab ls_pa_days_iv, mis
gen ls_pa_4days_iv=ls_pa_days_iv
recode ls_pa_4days_iv 7=0 6=0 5=0 4=0 3=1 2=1 1=1 0=1
lab var ls_pa_4days_iv "Imaging visit low PA - under 4 days/week moderate physical activity >10 mins"
lab val ls_pa_4days_iv ls_pa_4days_bllab

**vigorous physical activity- imaging visit
gen ls_vpa_days_iv=n_904_2_0 
recode ls_vpa_days_iv -7=. -1=. -3=.
lab var ls_vpa_days_iv "Imaging visit number of days/week vigourous physical activity >10 mins"
tab ls_vpa_days_iv, mis
gen ls_vpa_4days_iv=ls_vpa_days_iv
recode ls_vpa_4days_iv 7=0 6=0 5=0 4=0 3=1 2=1 1=1 0=1
lab var ls_vpa_4days_iv "Imaging visit low PA - under 4 days/week vigorous physical activity >10 mins"
lab val ls_vpa_4days_iv ls_vpa_4days_bllab




***alcohol use ***************************************

*Alcohol consumption frequency- baseline
gen ls_alc_freq_bl=n_1558_0_0
recode ls_alc_freq_bl -3=.
lab var ls_alc_freq_bl"Baseline alcohol intake frequency"
lab def ls_alc_freq_bllab 1"Daily or almost daily" 2"Three or four times a week" 3"Once or twice a week" 4"One to three times a month" 5"Special occasions only" 6"Never"
lab val ls_alc_freq_bl ls_alc_freq_bllab
gen ls_alc_freq3cats_bl=ls_alc_freq_bl
recode ls_alc_freq3cats_bl 5/6=0 3/4=1 1/2=2
lab var ls_alc_freq3cats_bl"Baseline alcohol intake frequency - 3 categories"
lab def ls_alc_freq3cats_bllab 0"Rarely or never" 1"1-8 times per month" 2"16 times per month- every day" 
lab val ls_alc_freq3cats_bl ls_alc_freq3cats_bllab
gen ls_bin_alc_freq_bl=ls_alc_freq_bl
replace ls_bin_alc_freq_bl=0 if ls_alc_freq_bl!=1
replace ls_bin_alc_freq_bl=. if ls_alc_freq_bl==.
tab ls_bin_alc_freq_bl
lab var ls_bin_alc_freq_bl"Alcohol intake daily/ almost vs. not"
lab def ls_bin_alc_freq_bllab 0"Alcohol intake less than daily" 1"Alcohol intake daily/ almost"
lab val ls_bin_alc_freq_bl ls_bin_alc_freq_bllab

*Alcohol consumption frequency- imaging visit
gen ls_alc_freq_iv=n_1558_2_0
recode ls_alc_freq_iv -3=.
lab var ls_alc_freq_iv"Imaging visit alcohol intake frequency"
lab def ls_alc_freq_ivlab 1"Daily or almost daily" 2"Three or four times a week" 3"Once or twice a week" 4"One to three times a month" 5"Special occasions only" 6"Never"
lab val ls_alc_freq_iv ls_alc_freq_bllab
gen ls_alc_freq3cats_iv=ls_alc_freq_iv
recode ls_alc_freq3cats_iv 5/6=0 3/4=1 1/2=2
lab var ls_alc_freq3cats_iv "Imaging visit alcohol intake frequency - 3 categories"
lab def ls_alc_freq3cats_ivlab 0"Rarely or never" 1"1-8 times per month" 2"16 times per month- every day" 
lab val ls_alc_freq3cats_iv ls_alc_freq3cats_ivlab
gen ls_bin_alc_freq_iv=ls_alc_freq_iv
replace ls_bin_alc_freq_iv=0 if ls_alc_freq_iv!=1
replace ls_bin_alc_freq_iv=. if ls_alc_freq_iv==.
tab ls_bin_alc_freq_iv
lab var ls_bin_alc_freq_iv"Alcohol intake daily/ almost vs. not"
lab def ls_bin_alc_freq_ivlab 0"Alcohol intake less than daily" 1"Alcohol intake daily/ almost"
lab val ls_bin_alc_freq_iv ls_bin_alc_freq_ivlab


*********smoking status*********************************************************************

**smoking status - baseline
gen smoking_status= n_20116_0_0 if n_20116_0_0>-1
label define smokelab 0"Never" 1"Former" 2"Current" 
label values smoking_status smokelab

gen ever_never=smoking_status
recode ever_never 2=1

**smoking status- imaging visit
gen smoking_status2= n_20116_2_0 if n_20116_2_0>-1
label values smoking_status2 smokelab





**genetic principal components

forvalues i = 1(1)20 {
rename n_22009_0_`i' pc`i'
}


****variable which indicates if participants have any MRI imaging data
gen any_MRI_data=1 if VAT!=. | ASAT!=. | ATMFI!=. | PTMFI!=. 



