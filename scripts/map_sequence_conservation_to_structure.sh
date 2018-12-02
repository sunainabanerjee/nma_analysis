#! /bin/bash


exedir="$(dirname $(readlink -f $0))"
consv_exe="${exedir}/conservation_from_sequence_alignment.sh"
restag_exe="${exedir}/map_seqpos2restag.sh"
poscorr_exe="${exedir}/map_struct2seqpos_corr.sh"



if [ ! -x "${consv_exe}" ] || [ ! -x "${restag_exe}" ] || [ ! -x "${poscorr_exe}" ]
then
  echo "Error: setup error!!" 1>&2; 
  exit 1
fi

usage(){ 
  echo "Usage: $0 -a <alignment file> -r <reference seq> 
                  -s <reference structure> [-S <reference chain>]  
                  -t <target structure> [-T <target chain>] [-c]" 1>&2; 
}

cleanup(){
  rm -f $tfile $rfile $qfile $xfile;
}


halt(){ 
  cleanup;
  exit 1; 
}


while getopts ":a:r:s:S:t:T:c" o;
do
  case "${o}" in
     r)
        tag="${OPTARG}"
        ;;
     a)
        alignment=$(readlink -f "${OPTARG}");
        ;;
     c)
        consv=1;
        ;;
     s)
        ref_struct=$(readlink -f "${OPTARG}");
        ;;
     S)
        ref_chain=${OPTARG};
        ;;
     t)
        tgt_struct=$(readlink -f "${OPTARG}");
        ;;
     T)
        ref_chain=${OPTARG};
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${alignment}" ] || [ ! -f "${alignment=}" ]
then
   usage;
   exit 1;
fi

if [ -z "${tag}" ] || [ $(grep -P "^>.*${tag}" "${alignment}"| wc -l) -ne 1 ]
then
  usage;
  exit 1;
fi

if [ -z "${ref_struct}" ] || [ -z "${tgt_struct}" ] || [ ! -f "${ref_struct}" ] || [ ! -f "${tgt_struct}" ]
then
   usage;
   exit 1;
fi

tfile=$(tempfile)
rfile=$(tempfile)
qfile=$(tempfile)
xfile=$(tempfile)

if [ ! -z ${consv} ]
then
  "${consv_exe}" -t "${tag}" -a "${alignment}" -c > $tfile

  if [ $? -ne 0 ];then  halt; fi
else
  "${consv_exe}" -t "${tag}" -a "${alignment}" > $tfile

  if [ $? -ne 0 ];then  halt; fi
fi
  

if [ ! -z ${ref_chain} ]
then
   "${restag_exe}" -p "${ref_struct}" -c "${ref_chain}" > $rfile
    if [ $? -ne 0 ];then  halt; fi
else
   "${restag_exe}" -p "${ref_struct}"  > $rfile
    if [ $? -ne 0 ];then  halt; fi
fi


if [ ! -z ${tgt_chain} ]
then
   "${restag_exe}" -p "${tgt_struct}" -c "${tgt_chain}" > $qfile
    if [ $? -ne 0 ];then  halt; fi
else
   "${restag_exe}" -p "${tgt_struct}"  > $qfile
    if [ $? -ne 0 ];then  halt; fi
fi

cmd="${poscorr_exe} -1 ${ref_struct} -2 ${tgt_struct} "

if [ ! -z  ${ref_chain} ]
then
   cmd="${cmd} -a $ref_chain";
fi

if [ ! -z ${tgt_chain} ]
then
   cmd="${cmd} -a $tgt_chain";
fi 

$cmd > $xfile
if [ $? -ne 0 ];then  halt; fi


awk -F, -vfile="${qfile}" -vtfile="${xfile}" 'BEGIN{ 
                              while( getline line < file){
                                split(line, arr, ",");
                                if( length(arr) == 2)
                                    map[arr[1]] = arr[2];
                              }
                              while( getline line < tfile ){
                                split(line, arr, ",");
                                if( length(arr) == 3)
                                   corr[arr[2]] = arr[3];
                              }
                            } NF==3 && corr[$2] != "" && map[corr[$2]] != "" {
                               printf("%s,%.2f\n",map[corr[$2]],$3);
                            }' "$tfile"

cleanup;
