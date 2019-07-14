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

./run_all_empPvals.sh # run the sample_meanCorr_empPvals.R


Rscript empPvalFC_empPvalCorr.R








