#!/usr/bin/bash

BAM=$1
REFDIR=$2
T1000=$3
SEX=$4
NORMAL_COV=$5
THREADS=$6
NNUM=$7

BIN="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export LD_LIBRARY_PATH=${BIN}/../Weaver_SV/lib:${LD_LIBRARY_PATH}

if [[ "${NNUM}" -eq 1 ]]; then
    SUFF=_num
fi

if ! [[ ${THREADS} =~ $re ]] ; then
    THREADS=8
fi

GAP=${BIN}/../data/GAP_20140416_num
GAPALPHA=${BIN}/../data/GAP_20140416

ARR=(WIG SV SNP)
for K in ${ARR[@]}; do
    echo $K,${ARR}
    perl $BIN/Weaver_pipeline.pl ALL ${K} \
        -f ${REFDIR} \
        -g ${GAPALPHA} -b ${BAM} \
        -k ${T1000} -s ${SEX} \
        -p $THREADS >weaver_${K} 2>weaver_${K}_error&
done
wait

FASTA=${REFDIR}.fasta
if [ ! -f ${FASTA} ]; then
    FASTA=${REFDIR}.fa
fi

${BIN}/Weaver PLOIDY \
    -f ${FASTA} \
    -s SNP_dens \
    -S ${BAM}.Weaver.GOOD -g ${GAP%$SUFF} \
    -m ${BIN}/../data/wgEncodeCrgMapabilityAlign100mer_number.bd \
    -w ${BAM}.wig -r 1 \
    -p ${THREADS} >weaver_ploidy 2>weaver_ploidy_error
    #-t ${CANCER_COV} -n ${NORMAL_COV} \

$BIN/solo_ploidy TARGET 2 > coverage_estimates        

CANCER_COV=`grep "Tumor Haplotype Coverage" coverage_estimates| cut -d ':' -f 2`  #`cat ${BAM}.coverage`
NORMAL_COV=`grep "Normal Haplotype Coverage" coverage_estimates | cut -d ':' -f 2`  #`cat ${BAM}.coverage`

re='^[0-9]+([.][0-9]+)?$'
#if ! [[ ${NORMAL_COV} =~ $re ]] ; then
#    NORMAL_COV=0
#fi
#NORMAL_COV=0

${BIN}/Weaver LITE \
    -f ${FASTA} \
    -s SNP_dens \
    -S ${BAM}.Weaver.GOOD -g ${GAP%$SUFF} \
    -m ${BIN}/../data/wgEncodeCrgMapabilityAlign100mer_number.bd \
    -w ${BAM}.wig -r 1 \
    -t ${CANCER_COV} -n ${NORMAL_COV} \
    -p  ${THREADS} >weaver_lite 2>weaver_lite_error
