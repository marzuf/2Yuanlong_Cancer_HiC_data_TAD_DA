#!/usr/bin/bash

# ./all_cmds.sh

Rscript create_sample_around_TADs_fixKb_twosided.R

Rscript create_sample_around_TADs_sameKb_twosided.R

Rscript create_sample_around_TADs_sameNbr_twosided.R


Rscript coexpr_around_TADs.R fixKb 1000000

Rscript coexpr_around_TADs.R sameKb

Rscript coexpr_around_TADs.R sameNbr


Rscript look_sampling_features.R

Rscript look_sampling_corr.R

Rscript check_left_and_right.R

#./run_all_empPvals.sh # run the sample_meanCorr_empPvals.R # modified, can be directly run as:
Rscript sample_meanCorr_empPvals.R


Rscript empPvalFC_empPvalCorr.R

# > run step19 pipeline in Yuanlong_Cancer_HiC_data_TAD_DA

Rscript SAM_emp_measurement_meanCorr.R

Rscript SAM_emp_measurement_meanCorr_allDS.R

Rscript scatterplot_FDR_byDS.R
Rscript scatterplot_FDR_allDS.R

Rscript meanCorr_meanFC_byCmpType.R

Rscript meanCorr_meanFC_FDR.R

Rscript boxplot_FDRthresh_nbrSignif.R


Rscript combined_pvals.R


Rscript combined_pvals_allDS.R


Rscript create_empPvalComb.R


Rscript intersect_FDR_combPval.R


Rscript features_signif_TADs.R

Rscript plot_signif_TADs.R

Rscript intersect_FDR_combPval_lolliPlot.R


Rscript cmp_same_conditions.R
Rscript geneRank_vs_tadRank.R

Rscript boundary_effect_gene.R

Rscript boundary_effect_TAD.R

Rscript combined_pvals_meth_cmp.R

Rscript create_sample_alternative.R
Rscript coexpr_alternative_samplings.R
Rscript cmp_coexpr_dist.R





