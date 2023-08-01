#!/bin/zsh

# !DEPRECATED

SCRIPT_NAME=$0

helpFunction()
{
   echo ""
   echo "Usage: $SCRIPT_NAME -t bool"
   echo -e "\t-t bool (to remove tsv files)"
   exit 1 # Exit script after printing help
}

while getopts "t:h" opt
do
   case "$opt" in
      t ) rm_tsv="$OPTARG" ;;
      h ) helpFunction ;; 
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ $# -eq 0 ]; then
  helpFunction
fi

# Download all coffea data with valid GPS lat/lon from Madagascar in specific year ranges
# NO API key on GBIF since GitHub account and simple http auth on GBIF's API
# curl --i --user simlal:??? -X GET -o gbif_coff_mg_all.zip "https://api.gbif.org/v1/occurrence/download?country=MG&has_coordinate=true&has_geospatial_issue=false&taxon_key=2895315&occurrence_status=present"
# curl --i --user simlal:??? -X GET -o gbif_coff_mg_2000-2023.zip "https://api.gbif.org/v1/occurrence/download?country=MG&has_coordinate=true&has_geospatial_issue=false&taxon_key=2895315&year=2000,2023&occurrence_status=present"
# curl --i --user simlal:??? -X GET -o gbif_coff_mg_2010-2023.zip "https://api.gbif.org/v1/occurrence/download?country=MG&has_coordinate=true&has_geospatial_issue=false&taxon_key=2895315&year=2010,2023&occurrence_status=present"

# # # Move data to coffea dir
# cd ~/Downloads
# ls | grep gbif | xargs -I {} mv {} ../programming/python_projects/coffeaPhyloGeo/data/coffea_example/gbif_data

# Unzip and convert tsv to csv
cd ~/programming/python_projects/coffeaPhyloGeo/data/coffea_example/gbif_data/

zipfiles=(*.zip(N))
if [ ${#zipfiles} -gt 0 ]; then    # Safety check for .zipfiles in current dir
  for f in "${zipfiles[@]}"; do 
    tsv_file="${f%.zip}.tsv"
    csv_file="${f%.zip}.csv"
    unzip -p "$f" > $tsv_file;

    # Convert tsv-csv with fields enclosed in quotes
    awk 'BEGIN {OFS = FS = "\t"} {for(i=1; i<=NF; i++) $i="\""$i"\""; print $0}' "${tsv_file}" | tr '\t' ',' > "${csv_file}"
    rm "$f"

    if [ $rm_tsv = "true" ]; then
      rm "$tsv_file"
    fi
  
  done
else
  echo "No .zip files found in current directory."
fi