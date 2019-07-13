#!/usr/bin/bash

# ./run_all.sh

rscript="sample_meanCorr_empPvals.R"

all_ds=(

 # LIVER
"GSE105381_HepG2_40kb TCGAlihc_norm_lihc"
"GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1"

"GSE58752_liver_40kb TCGAlihc_norm_lihc"
"GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1"


#  BREAST

"ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"
"GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas"


#  KIDNEY
"ENCSR079VIJ_G401_40kb TCGAkich_norm_kich"
"ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"


#  SKIN
"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"

"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"

"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"

#  LUNG
"ENCSR444WCZ_A549_40kb TCGAluad_norm_luad"
"ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad"

"ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR"
"ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"

"ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker"
"ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker"

"ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS"
"ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS"

"ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc "
"ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc"

#  PANCREAS
"Panc1_rep12_40kb TCGApaad_wt_mutKRAS"


#  PROSTATE	
"ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad"
"GSE73782_PC3_40kb TCGAprad_norm_prad"
"GSE118514_RWPE1_40kb TCGAprad_norm_prad"

#  GBM
"GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal"
"GSE105194_cerebellum_40kb TCGAgbm_classical_neural"
"GSE105194_cerebellum_40kb TCGAgbm_classical_proneural"
"GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc"

"GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal"
"GSE105194_spinal_cord_40kb TCGAgbm_classical_neural"
"GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural"
"GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc"

 # COLORECTAL
"GSE105318_DLD1_40kb TCGAcoad_msi_mss"

 # LYMPHOBLAST
"K562_40kb TCGAlaml_wt_mutFLT3"

)


for ds in "${all_ds[@]}"; do
    
    echo Rscript $rscript $ds
    Rscript $rscript $ds
done




