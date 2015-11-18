#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Call discordance matrix parser on all available matrices
# ----------------------------------------------------------------------------------------

for disc in results/*discordance_matrix; do

    files=`echo $disc | sed -e "s:results/diff_::" -e "s/\.diff\.discordance_matrix//"`

    OIFS=$IFS
    IFS='_'
    read -ra file_arr <<< "$files"
    IFS=$OIFS
    
    file1=${file_arr[0]}
    file2=${file_arr[1]}

    echo "Processing comparison between $file1 and $file2..." 1>&2
    
    parsed=`Rscript scripts/parse_discordance_matrix.R $disc`
    
    echo $file1 $file2 $parsed
done
