#!/bin/bash

# Should be run in wbbi201 to use 80 cores
# Use screen if using putty
## run 5_prepare_and_run_main.sh in source parent folder were test was run
# use dos2unix

conda deactivate
echo "Going to source parent folder were test was run (in test folder)"
cd /artemis/com/shg29ny/2020_04_21_afum_ps
mkdir main

echo coping required files from test to main
cd main
cp ../test/*.sh .
cp -r ../test/github_uploaded_scripts .
cp -r ../test/Input_folder .

## The main folder should look like
## 1_collect_data_and_preprocess.sh  2_set_up_project_folders.sh  3_process_clusters_for_PAML_input.sh  4_cleanup_for_rerun_if_needed.sh  github_uploaded_scripts  Input_folder

echo Editing scripts - replacing test to main using sed to change the folder path
sed -i 's/test/main/g' *.sh
sed -i 's/maining/testing/g' 1_collect_data_and_preprocess.sh

echo running first script
source 1_collect_data_and_preprocess.sh &>> log_script1

echo running second script
source 2_set_up_project_folders.sh &>> log_script2

echo using sed to replace the input subset file parsed_clusters_tail with final file parsed_clusters_Final
sed -i 's/parsed_clusters_tail/parsed_clusters_Final/g' 3_process_clusters_for_PAML_input.sh

echo running third-main script
source 3_process_clusters_for_PAML_input.sh &>> log_script3









