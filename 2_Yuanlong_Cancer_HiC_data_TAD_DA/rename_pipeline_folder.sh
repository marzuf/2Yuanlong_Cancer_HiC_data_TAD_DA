#!/usr/bin/bash

# ./rename_pipeline_folder.sh

pipFolder="../Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER"


#oldName="14f210000_cumulAllDown_limited_AUC"
#newName="14f210000_cumulAllDown_limited_AUC_RENAMED"


#oldName="19_SAM_emp_measurement"
#newName="1910000_SAM_emp_measurement"

oldName="11sameNbr_runEmpPvalCombined"
newName="11sameNbr10000_runEmpPvalCombined"


all_folders=( $(realpath $pipFolder/*/*/$oldName) )

for folder in ${all_folders[@]}; do

    echo "-> mv $folder `dirname $folder`/$newName"
    mv $folder `dirname $folder`/$newName

done 
