cmp_names <- c(

    "wt_vs_mut" = "wt. vs. mut.",
    "norm_vs_tumor" = "normal vs. tumor",
    "subtypes" = "subtypes"

) 

hicds_names <- c( 
"GM12878_40kb" ="GM12878")         

exprds_names <- c(
  "TCGAbrca_lum_bas" = "BRCA luminal vs. basal",             
  "TCGAkich_norm_kich" = "KICH normal vs. tumor",           
 "TCGAskcm_lowInf_highInf" = "SKCM low inf. vs. high inf.", 
 "TCGAskcm_wt_mutBRAF" = "SKCM BRAFwt vs. BRAFmut",          
 "TCGAskcm_wt_mutCTNNB1" = "SKCM CTNNB1wt vs. CTNNB1mut",          
 "TCGAprad_norm_prad" = "PRAD normal vs. tumor",           
 "TCGAluad_mutKRAS_mutEGFR" = "LUAD KRASmut vs. EGFRmut",
 "TCGAluad_nonsmoker_smoker" = "LUAD nonsmoker vs. smoker",
 "TCGAluad_norm_luad" = "LUAD normal vs. tumor",
 "TCGAluad_wt_mutKRAS" = "LUAD KRASwt vs. KRASmut",          
 "TCGAlusc_norm_lusc" = "LUSC normal vs. tumor",
 "TCGAcoad_msi_mss" = "COAD MSI vs. MSS",  
 "TCGAgbm_classical_mesenchymal" = "GBM classical vs. mesenchymal",
 "TCGAgbm_classical_neural"      = "GBM classical vs. neural",
 "TCGAgbm_classical_proneural"  = "GBM classical vs. proneural",
 "TCGAlgg_IDHwt_IDHmutnc" = "LGG IDHwt vs. IDHmut",
 "TCGAlihc_norm_lihc" = "LIHC normal vs. tumor",
 "TCGAlihc_wt_mutCTNNB1" = "LIHC CTNNB1wt vs. CTNNB1mut",
 "TCGApaad_wt_mutKRAS" = "PAAD KRASwt vs. KRASmut",
 "TCGAlaml_wt_mutFLT3" = "LAML FLT3wt vs. FLT3mut"
)

cond1_names <- sapply(exprds_names, function(x) gsub(".+ (.+) vs. (.+)", "\\1", x))
names(cond1_names) <- sapply(names(exprds_names), function(x) gsub("TCGA.+_(.+)_(.+)", "\\1", x))
cond2_names <- sapply(exprds_names, function(x) gsub(".+ (.+) vs. (.+)", "\\2", x))
names(cond2_names) <- sapply(names(exprds_names), function(x) gsub("TCGA.+_(.+)_(.+)", "\\2", x))



