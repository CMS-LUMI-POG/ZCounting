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

# need to shift arguments to not confuse other scripts
shift
shift
shift
shift

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


WEBBASE=/eos/home-c/cmszcont/www/private/Run3/                # CMS ZCounting service account (cmszcont)
#WEBBASE=/eos/user/t/tmenezes/www/private/Run3/test_runs        # Your personal account
WEBSUB=$YEAR/
WEBDIR=$WEBBASE/$WEBSUB/$normtag/$jsonType/$etaRegion  # Include normtag in WEBDIR path

mkdir -p $WEBBASE/$WEBSUB
mkdir -p $WEBBASE/$WEBSUB/$normtag
mkdir -p $WEBBASE/$WEBSUB/$normtag/$jsonType
mkdir -p $WEBBASE/$WEBSUB/$normtag/$jsonType/$etaRegion
mkdir -p $WEBDIR  # Create the output directory for the normtag


CMSSWBASE=/afs/cern.ch/user/c/cmszcont/CMSSW_12_4_18/src/                 # CMS ZCounting service account (cmszcont)
#CMSSWBASE=/afs/cern.ch/user/t/tmenezes/work/private/CMSSW_12_4_18/src/     # Your personal account
ZHARVESTER=$CMSSWBASE/ZCounting/ZHarvester

PDs=('/SingleMuon/' '/Muon/' '/Muon0/' '/Muon1/')  # look for all possible primary datasets


# During the data-taking the new data will be stored in the DQMOffline
#COREDIR="https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run${YEAR}/"


# After the data taking, the data will be stored in DQMGUI_data
COREDIR="/eos/cms/store/group/comm_dqm/DQMGUI_data/Run${YEAR}"
#COREDIR="/eos/user/t/tmenezes/ZCounting_Data/DQM_Histograms/store/group/comm_dqm/DQMGUI_data/Run${YEAR}"


PDDIR="$COREDIR$PDDIR/"




RESOURCES_OUTPUT="/eos/home-c/cmszcont/data/zcounting"          # CMS ZCounting service account (cmszcont)
#RESOURCES_OUTPUT="/eos/user/t/tmenezes/ZCounting_Data/"          # Your personal account
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
#voms-proxy-init --voms cms -rfc
CERT=`voms-proxy-info --path`
echo $CERT
#CERT=`voms-proxy-info --path`

##########################################################
## download normtag and golden/DCS json, produce byls csv
##########################################################

echo 'Set brilcalc environment'
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
echo "BRILCALC_OUTPUT is set to: ${BRILCALC_OUTPUT}"

# define the normtag path
normtagPath="/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_${normtag}.json"

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

	# Path using the Normtags from the normtag directory
        #brilcalcFileName="${lumimaskFileName/.json/.csv}"

	# Writing the normtag directly (i.e. hfoc23v07)
	brilcalcFileName="${lumimaskFileName/.json/}_${normtag}.csv"
        brilcalcFilePath=${BRILCALC_OUTPUT}/${brilcalcFileName}


        if [[ -e $brilcalcFilePath ]]; then
            echo "File ${brilcalcFilePath} exists already"
        else
            echo "Make brilcalc from ${lumimaskName} to ${brilcalcFilePath}"
            cp $lumimaskFilePath ./$lumimaskName

        # definition of brilcalc alias found via `alias brilcalc`:
        if [[ "$jsonType" == "DCSOnly_TkPx" ]]; then
            singularity -s exec  --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest \
            brilcalc lumi --datatag online -b "STABLE BEAMS" --byls -u /fb -i $lumimaskName -o $brilcalcFileName
        else
            singularity -s exec --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest \
	    brilcalc lumi --normtag $normtag -b "STABLE BEAMS" --byls -u /fb -i $lumimaskName -o $brilcalcFileName
        fi

	# brilcalc lumi --normtag $normtagPath -b "STABLE BEAMS" --byls -u /fb -i $lumimaskName -o $brilcalcFileName
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
rm -rf $filelist


echo "COREDIR: $COREDIR"


for index in "${!PDs[@]}"
do
    PD=${PDs[$index]}
    PDDIR="$COREDIR$PD"
    echo "PDDIR: $PDDIR"

    echo "Retrieve contents from $PDDIR"
    # Check if the directory exists
    if [ -d "$PDDIR" ]; then
        # Recursively list files in the directory and append only .root files to the filelist log
        #find "$PDDIR" -type f -name "*.root" >> "$filelist"
	    find "$PDDIR" -type f -name "*.root" | grep "PromptReco-v1" >> "$filelist"
    else
        echo "Warning: Directory $PDDIR does not exist" >&2
    fi
done

echo "File listing completed. Check the log file at: $filelist"


######################################################
# Skip copying files and directly use the file paths
######################################################

#echo "Using files in $filelist"
#ls -l $filelist


echo "Copy files in $filelist"
ls -l $filelist
cat $filelist | while read LINE

do
    # Determine the relative path from the source file
    RELATIVE_PATH=${LINE#/eos/cms/}

    # Construct the target file path, preserving the folder structure
    TARGET=${HISTOGRAM_OUTPUT}/${RELATIVE_PATH}
    #TARGET=${HISTOGRAM_OUTPUT}

    # Ensure the target directory exists
    TARGET_DIR=$(dirname $TARGET)
    mkdir -p $TARGET_DIR

    # Check if the file already exists in the output directory
    if [ -e $TARGET ]; then
        echo "File exists already: $TARGET"
    else
        echo "Copying file: $LINE to $TARGET"
        #cp $LINE $TARGET  # Use appropriate copy command (e.g., cp, xrdcp, eos cp)
    fi
done

#############################################
## Run ZCounting fits
#############################################

echo "Fit files in $filelist"
cat $filelist | while read LINE

do
    path="$LINE"  # Directly use the line read from the filelist
    if ! [[ -e $path ]]; then
        echo "ERROR: file $path does not exist but it should"
    fi
    ls $path
    echo now at path $path


    # extract the run number from the filename
    old_IFS=$IFS # Save the current value of IFS, set the field separator to "_"
    IFS="_"
    filename=$(basename "$path")
    read -ra parts <<< "$filename" # Convert the string to an array
    IFS=$old_IFS # Restore the value of IFS to its original value
    run_number="${parts[2]:4}" # Access the third item (index 2) from the array


    # Update path to the selected file and ensure only this file is processed
    #path=$selected_file
    echo "Selected file: $path"

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
    echo input path ${path}
    echo histogram output ${HISTOGRAM_OUTPUT}
    echo brilcalc $BRILCALC_OUTPUT
    echo output $FITS_OUTPUT
    ./ZCounting_DQM.py -b $beginRun -e $endRun -i ${COREDIR} \
    --byLsCSV $byLsCSV -o $FITS_OUTPUT --etaMin $etaMin --etaCut $etaMax 
done




echo "Produce plots per fill"
#python3 Plotting/plot_ZLumi_RefLumi.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Fills --year ${YEAR} --xsec ${xsec} --rrange 0.9 1.1

echo "Produce stability plots"
#python3 Plotting/plot_stability.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Stability --year ${YEAR}

echo "Produce linearity plots"
#echo csvfile: $CSVS_OUTPUT/Mergedcsvfile_perMeasurement.csv webdir: $WEBDIR and xsec: $xsec 
#python3 Plotting/plot_linearity.py -r ${CSVS_OUTPUT}/Mergedcsvfđxsecile_perMeasurement.csv -o ${WEBDIR}/Linearity --year ${YEAR} --xsec ${xsec}

# TODO: plot_ZSummary.py
echo "Produce Z Summary plots"
#python3 Plotting/plot_ZSummary.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Summary --year ${YEAR} --xsec ${xsec}
python3 Plotting/plot_ZSummary.py -r ${CSVS_OUTPUT}/Mergedcsvfile_perMeasurement.csv -o ${WEBDIR}/Summary --year ${YEAR} --xsec ${xsec} --saveas ${YEAR}_zcount

#echo Running ZCounting_DQM for $beginRun and $endRun with the input from $COREDIR the byLsCSV $byLsCSV with the output in $FITS_OUTPUT
