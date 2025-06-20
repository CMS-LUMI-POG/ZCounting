#!/bin/bash
echo "Awake"
. ~/.bashrc

# input settings
if [[ $# -ne 4 ]]; then
    echo "Usage: get_brilcalc <year> <json type> <eta region> <normtag>"
    return 1
fi

YEAR=$2       # can be 2023, 2024, ...
jsonType=$1   # can be Golden or DCSOnly_TkPx
etaRegion=$3  # can be Barrel or Inclusive
normtag=$4    # can be BRIL, HFET, etc.

normtag="${normtag//$'\r'/}"
normtag="${normtag//[[:space:]]/}"

# need to shift arguments to not confuse other scripts
shift
shift
shift
shift

allowed_normtags=("hfet" "dt" "hfoc" "pxl" "pltzero" "bcm1futca" "bcm1f")

if [[ $etaRegion == "Barrel" ]]; then 
    etaMin=0.0
    etaMax=0.9
    xsec=170 # for 13.6 TeV (from 2022 E,F,G)
elif [[ $etaRegion == "Endcap" ]]; then 
    etaMin=0.9
    etaMax=2.4
else
    etaMin=0.0
    etaMax=2.4
    xsec=682 # for 13.6 TeV (from 2022 E,F,G)
fi

echo "Processing data from $YEAR using ${jsonType}.json files for ${etaRegion} detector region and normtag ${normtag}"

DATA=/eos/home-c/cmszcont/data

WEBBASE=/eos/home-c/cmszcont/www/private/Run3/
WEBSUB=$YEAR/
WEBDIR=$WEBBASE/$WEBSUB/$normtag/$jsonType/$etaRegion  # Include normtag in WEBDIR path

mkdir -p $WEBBASE/$WEBSUB
mkdir -p $WEBBASE/$WEBSUB/$normtag
mkdir -p $WEBBASE/$WEBSUB/$normtag/$jsonType
mkdir -p $WEBBASE/$WEBSUB/$normtag/$jsonType/$etaRegion
mkdir -p $WEBDIR  # Create the output directory for the normtag

CMSSWBASE=/afs/cern.ch/user/c/cmszcont/CMSSW_12_4_18/src/
ZHARVESTER=$CMSSWBASE/ZCounting/ZHarvester

PDs=('/SingleMuon/' '/Muon/' '/Muon0/' '/Muon1/')  # look for all possible primary datasets

COREDIR="https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run${YEAR}/"
PDDIR="$COREDIR$PDDIR/"

RESOURCES_OUTPUT="/eos/home-c/cmszcont/data/zcounting/"
HISTOGRAM_OUTPUT="${RESOURCES_OUTPUT}/DQM_Histograms/"
LUMIMASK_OUTPUT="${RESOURCES_OUTPUT}/Lumimasks/"
BRILCALC_OUTPUT="${RESOURCES_OUTPUT}/Brilcalc_byLsCSV/${YEAR}/${normtag}/${jsonType}"

FITS_OUTPUT="${WEBDIR}/RunsAndFits/"
CSVS_OUTPUT="$FITS_OUTPUT/csvFiles"

mkdir -p $FITS_OUTPUT
mkdir -p $CSVS_OUTPUT
mkdir -p $BRILCALC_OUTPUT

echo "Setup environment"

# get grid certificate
voms-proxy-init --voms cms -rfc
CERT=`voms-proxy-info --path`
echo $CERT

##########################################################
## download normtag and golden/DCS json, produce byls csv
##########################################################

echo 'Set brilcalc environment'
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env

# Define how to pass the normtag: via --type (for detector names) or --normtag (for specific normtag)
#if [[ " ${allowed_normtags[@]} " =~ " ${normtag} " ]]; then
if [[ "${normtag^^}" == "BRIL" ]]; then
    echo "Using BRIL online lumi"
    normtagOption="--datatag online"
elif [[ " ${allowed_normtags[*]} " == *" ${normtag} "* ]]; then
    echo "Using detector type: ${normtag}"
    normtagOption="--type ${normtag}"
else
    echo "Using specific normtag argument: ${normtag}"
    normtagOption="--normtag ${normtag}"
fi


lumimaskbase="cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions${YEAR:2:4}/"
if [[ "$jsonType" == "DCSOnly_TkPx" ]]; then
    lumimaskbase=$lumimaskbase/DCSOnly_JSONS/dailyDCSOnlyJSON/
fi

echo "Get all ${jsonType} jsons"
lumimask="https://${lumimaskbase}"

curl -k --cert $CERT -X GET $lumimask | awk -F "<a href=\"" '{ print $2} ' | while read -r piece; do
    echo "Now at piece $piece"
    if [[ "${piece,,}" == *"${jsonType,,}.json"* ]]; then
        modified_piece=${piece%%</a>*}   # Removes </a>
        modified_piece=${modified_piece%%\">*}   # Removes ">
        lumimaskFileName=${modified_piece}
        lumimaskName=${lumimaskFileName}
        lumimaskFilePath=$LUMIMASK_OUTPUT/${lumimaskbase}${lumimaskFileName}

        if [[ -e $lumimaskFilePath ]]; then
            echo "File ${lumimaskFilePath} exists already"
        else
            echo "Download ${lumimask}/${lumimaskFileName}"
            wget -r -k -4 ${lumimask}/${lumimaskFileName} -P $LUMIMASK_OUTPUT --certificate=$CERT --ca-certificate=$CERT 
        fi

        brilcalcFileName="${lumimaskFileName/.json/.csv}"
        brilcalcFilePath=${BRILCALC_OUTPUT}/${brilcalcFileName}


        if [[ -e $brilcalcFilePath ]]; then
            echo "File ${brilcalcFilePath} exists already"
        else
            echo "Make brilcalc from ${lumimaskName} to ${brilcalcFilePath}"
            cp $lumimaskFilePath ./$lumimaskName

        #definition of brilcalc alias found via `alias brilcalc`: 
        #if [[ "$jsonType" == "DCSOnly_TkPx" ]]; then
        #    singularity -s exec --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest brilcalc lumi --datatag online -b "STABLE BEAMS" --byls -u /fb -i $lumimaskName -o $brilcalcFileName
        #else
        singularity -s exec --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest brilcalc lumi $normtagOption -b "STABLE BEAMS" --byls -u /fb -i "$lumimaskName" -o "$brilcalcFileName"
        #fi

	    mv $brilcalcFileName $brilcalcFilePath
	    rm $lumimaskName
	    echo "brilcalc file created: ${brilcalcFileName}"
        fi
    fi
done


####################
# Get List Of Files
####################

cd $CMSSWBASE
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd $ZHARVESTER

echo "FITS_OUTPUT is set to: ${FITS_OUTPUT}"

# Create the log files with the files from the PDDIR (uses grid certs. and grid passwd)
echo "Get List of files"
filelist="${FITS_OUTPUT}/filelist.log"
#rm -rf $filelist

for index in "${!PDs[@]}"
do
    PD=${PDs[$index]}
    PDDIR=$COREDIR$PD

    echo "Attempting to retrieve content from $PDDIR"    
    command="curl -L -k --cert $CERT -X GET $PDDIR"
    echo $command
    echo "Output from curl for PDDIR: $(curl -L -k --cert $CERT -X GET $PDDIR)"
    for rundir in `$command | awk 's=index($0,"/000") { print substr($0,s+1,10)}'`
    do
        echo "Now at run $rundir"
        # loop over all files of each directory and write the file path in the filelist
        echo "$(curl -L -k --cert $CERT -X GET $PDDIR$rundir)" | awk -F '<tr><td>' '{ print  substr($2,60) } ' | while read -r piece; do
            if [ -n "$piece" ]; then
                echo $COREDIR${piece%%\'*} >> $filelist
            fi
        done
    done
done


######################################################
# Check which files exist already and download others
######################################################

echo "Copy files in $filelist"
ls -l $filelist
cat $filelist | while read LINE
do      
    if ! [ -e ${HISTOGRAM_OUTPUT}${LINE:8:1000} ]; then
        wget -r -k -4 $LINE -P $HISTOGRAM_OUTPUT --certificate=$CERT --ca-certificate=$CERT 
        # curl -k --cert $LINE -o $HISTOGRAM_OUTPUT
    else
        echo "File exists already"
        echo $LINE
    fi
done


#############################################
## Run ZCounting fits
#############################################


echo "Fit files in $filelist"
cat $filelist | while read LINE
do
    path=${HISTOGRAM_OUTPUT}${LINE:8:1000}
    if ! [[ -e $path ]]; then
        echo "ERROR: file $path does not exist but it should"
    fi

    echo now at path $path

    if [[ $path == *"_Muon1_"* ]]; then
        # loop only over Muon0, skip Muon1 but check if Muon1 exists
        continue
    elif [[ $path == *"_Muon0_"* ]]; then
        path1="${path//"Muon0"/"Muon1"}"
        if ! [ -e "$path1" ]; then
            echo "File ${path1} does not exist"
            continue
        fi
    fi

    # extract the run number from the filename
    old_IFS=$IFS # Save the current value of IFS, set the field separator to "_"
    IFS="_"
    filename=$(basename "$path")
    read -ra parts <<< "$filename" # Convert the string to an array
    IFS=$old_IFS # Restore the value of IFS to its original value
    run_number="${parts[2]:4}" # Access the third item (index 2) from the array

    # skip runs that are already processed
    if [[ -e "${CSVS_OUTPUT}/csvfile${run_number}.csv" ]]; then
        echo "Result csv file for run ${run_number} exists already, continue with next one"
        continue
    fi

    # skip some run ranges to save time
    if (( run_number > 347687 && run_number <= 355065 )) \
        || (( run_number > 362180 && run_number <= 362349 )) \
        || (( run_number > 362760 && run_number <= 366364 )) \
        || (( run_number > 372415 && run_number <= 373076 ))
    then
        continue
    fi

    # Define variables
    byLsCSV=$BRILCALC_OUTPUT
    beginRun=$(($run_number))
    endRun=$(($run_number + 1))

    # Get list of files in the directory
    files=$(ls "$BRILCALC_OUTPUT")

    # Filter out files containing "_era"
    filtered_files=()
    for file in $files; do
        if [[ ! "$file" == *"_era"* ]]; then
            filtered_files+=("$file")
        fi
    done

    # Filter files based on beginRun and endRun
    byLsCSVs=()
    for file in "${filtered_files[@]}"; do
        unset numbers
        parts=($(echo "$file" | tr "_" "\n"))
        for part in "${parts[@]}"; do
            if [[ "$part" =~ ^[0-9]+$ ]]; then
                numbers+=("$part")
            fi
        done
        if [[ "${#numbers[@]}" -ge 2 ]]; then
            start_num="${numbers[0]}"
            end_num="${numbers[1]}"
            if (( beginRun >= start_num && endRun <= end_num )); then
                byLsCSVs+=("$file")
            fi
        fi
        unset numbers
    done

    # Check if any files match and log the appropriate message or raise an error
    if (( ${#byLsCSVs[@]} > 0 )); then
        echo "byLS file ${byLsCSVs[0]} found"
        byLsCSV="$byLsCSV/${byLsCSVs[0]}"
    else
        echo "No byLS file found for run $run_number, exit"
        continue
    fi

    echo now at Run $run_number
    echo input ${HISTOGRAM_OUTPUT}${COREDIR:8}
    echo brilcalc $BRILCALC_OUTPUT
    echo output $FITS_OUTPUT
    ./ZCounting_DQM.py -b $beginRun -e $endRun -i ${HISTOGRAM_OUTPUT}${COREDIR:8} \
    --byLsCSV $BRILCALC_OUTPUT -o $FITS_OUTPUT --etaMin $etaMin --etaCut $etaMax 
done

#Fixing Mergedcsvfile_perMeasurement.csv creation
# ———  ———
CSVS_OUTPUT="$FITS_OUTPUT/csvFiles"
out="$CSVS_OUTPUT/Mergedcsvfile_perMeasurement.csv"
shopt -s nullglob
files=($CSVS_OUTPUT/csvfile*.csv)
if (( ${#files[@]} > 0 )); then
  echo "Merging ${#files[@]} per‐run CSVs into $out"
  head -n1 "${files[0]}" > "$out"
  for f in "${files[@]}"; do
    tail -n+2 "$f" >> "$out"
  done
else
  echo "No per‐run CSVs found; merged CSV stays empty"
fi
shopt -u nullglob
# ————————————————————————————————————————————————


#############################################
## Run ZCounting fits
#############################################

echo "Produce plots per fill"
python3 Plotting/plot_ZLumi_RefLumi.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Fills --year ${YEAR} --xsec ${xsec} --rrange 0.9 1.1 --normtag ${normtag}

echo "Produce stability plots"
python3 Plotting/plot_stability.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Stability --year ${YEAR}

# echo "Produce linearity plots"
# python3 Plotting/plot_linearity.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Linearity

# TODO: plot_ZSummary.py

echo "Produce stability plots"
python3 Plotting/plot_ZSummary.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Summary --year ${YEAR} --xsec ${xsec} --normtag ${normtag}


export PATH="${PATH}:${HOME}/php-plots/bin"
pb_copy_index.py ${WEBDIR} --recursive