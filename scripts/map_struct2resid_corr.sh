#!/bin/bash


exedir=$(dirname $(readlink -f $0))
struct2seq="${exedir}/map_struct2seqpos_corr.sh"
seq2res="${exedir}/map_resid2seqpos.sh"

if [ ! -x "${struct2seq}" ] || [ ! -x "${seq2res}" ]
then
  echo "Error: Set up issue, missing dependent scripts" 1>&2;
  exit 1;
fi

usage(){
  echo "Usage: $0 -1 <pdb file1> -2 <pdb file2> -a <pdb file1 chain id> -b <pdb file-2 chain id>" 1>&2;
}

while getopts ":1:2:a:b:" o;
do
  case "${o}" in
    1) pdbfile1=$(readlink -f ${OPTARG});
       ;;
    2) pdbfile2=$(readlink -f ${OPTARG});
       ;;
    a) chain1=${OPTARG};
       ;;
    b) chain2=${OPTARG};
       ;;
    *) usage;
       exit 1;
  esac
done
shift $((OPTIND-1))

if [ -z "${pdbfile1}" ] || [ ! -f "${pdbfile1}" ] || [ $(awk '/^ATOM/' "${pdbfile1}"| wc -l) -lt 4 ]
then
  usage;
  exit 1
fi

if [ -z "${pdbfile2}" ] || [ ! -f "${pdbfile2}" ] || [ $(awk '/^ATOM/' "${pdbfile2}"| wc -l) -lt 4 ]
then
  usage;
  exit 1
fi

if [ -z "${chain1}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile1}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain1=$(find_neighbor_pdb -f "${pdbfile1}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi

if [ -z "${chain2}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile2}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain2=$(find_neighbor_pdb -f "${pdbfile2}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi

if [ $(echo -n "$chain1" | wc -c) -ne 1 ] || [ $(echo -n "$chain2" | wc -c) -ne 1 ]
then
  usage;
  exit 1;
fi

tfile1=$(tempfile)
tfile2=$(tempfile)
cfile=$(tempfile)


${seq2res} -p "${pdbfile1}" -c ${chain1} > "${tfile1}"
${seq2res} -p "${pdbfile2}" -c ${chain2} > "${tfile2}"
${struct2seq} -1 "${pdbfile1}" -2 "${pdbfile2}" -a $chain1 -b $chain2 > "${cfile}"

awk -F, -vfile1=$tfile1 -vfile2=$tfile2 'BEGIN{  
                                           while( getline line < file1 ){
                                             split(line, arr, ",");
                                             if(length(arr) == 2)
                                                map1[arr[1]] = arr[2];
                                           }
                                           while( getline line < file2 ){
                                             split(line, arr, ",");
                                             if(length(arr) == 2)
                                                map2[arr[1]] = arr[2];
                                           }
                                         }
                                         NF==3{
                                           printf("%d,%d,%d\n",$1,map1[$2],map2[$3]);
                                         }' $cfile

rm -f $tfile1 $tfile2 $cfile

