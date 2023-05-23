#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

ZSTD=""
REGIONS=""
TARGETS=""
SAMPLES=""
ZSTD_LEVEL=""
unset -v NO_KEEP

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -f|--filename)
    FILENAME="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--regions)
    REGIONS="-r $2"
    shift # past argument
    shift # past value
    ;;
    -t|--targets)
    TARGETS="-t $2"
    shift # past argument
    shift # past value
    ;;
    -s|--samples)
    SAMPLES="-s $2"
    shift # past argument
    shift # past value
    ;;
    --zstd)
    ZSTD="--zstd"
    shift # past argument
    ;;
    --zstd-level)
    ZSTD_LEVEL="--zstd-level $2"
    shift # past argument
    shift # past value
    ;;
    --no-keep)
    NO_KEEP="YES"
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${FILENAME}" ]
then
    echo "Specify a filename with --filename, -f <filename>"
    exit 1
fi

echo "FILENAME        = ${FILENAME}"

TMPDIR=$(mktemp -d -t xsi_XXXXXX) || { echo "Failed to create temporary directory"; exit 1; }

echo "Temporary director : ${TMPDIR}"
echo "Region : ${REGIONS}"
echo "Targets : ${TARGETS}"
echo "Samples : ${SAMPLES}"

function exit_fail_rm_tmp {
    echo "Removing directory : ${TMPDIR}"
    rm -r ${TMPDIR}
    exit 1
}

# --variant-block-length 65536
# --variant-block-length 1024
"${SCRIPTPATH}"/../../xsqueezeit -c ${ZSTD} ${ZSTD_LEVEL} --maf 0.002 -f ${FILENAME} -o ${TMPDIR}/compressed.xsi || { echo "Failed to compress ${FILENAME}"; exit_fail_rm_tmp; }
"${SCRIPTPATH}"/../../xsqueezeit -x ${REGIONS} ${TARGETS} ${SAMPLES} -f ${TMPDIR}/compressed.xsi -o ${TMPDIR}/uncompressed.bcf || { echo "Failed to uncompress ${FILENAME}"; exit_fail_rm_tmp; }

command -v bcftools || { echo "Failed to find bcftools, is it installed ?"; exit_fail_rm_tmp; }

echo
#echo "Diff between original file and uncompressed compressed file :"
#echo "We expect at least one line that differs, which is the bcftools_viewCommand="
#echo

# Diffing two large VCF/BCFs will take huge amount of memory since diff is not
# a streaming program, it will load everything in memory first...
# E.g., 1KGP3 chr20 -> about 9 GB (uncompressed view output) times 2 (two files)
#diff <(bcftools view ${FILENAME}) <(bcftools view ${TMPDIR}/uncompressed.bcf) | tee ${TMPDIR}/difflog.txt
diff <(bcftools view ${REGIONS} ${TARGETS} ${SAMPLES} ${FILENAME}) <(bcftools view ${TMPDIR}/uncompressed.bcf) > ${TMPDIR}/difflog.txt
DIFFLINES=$(wc -l ${TMPDIR}/difflog.txt | awk '{print $1}')
#echo $DIFFLINES
if [ ${DIFFLINES} -gt 4 ]
then
    if [ -z "${NO_KEEP}" ]
    then
        echo
        echo "[KO] The files differ, check out ${TMPDIR}/difflog.txt"
        exit 1
    else
        exit_fail_rm_tmp # For unit testing
    fi
    exit 1
else
    echo
    echo "[OK] The files are the same"
fi

rm -r $TMPDIR
exit 0

# sdiff <(./console_app -x -f temp_test/chr20_mini.bin -s "^NA12878") <(bcftools view ../Data/pbwt/chr20_mini.bcf -s "^NA12878")