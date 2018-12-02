#! /bin/bash

exedir=$(dirname $(readlink -f $0))

burstExe="${exedir}/burst_model.sh"
sasaExe="${exedir}/calculate_sasa.sh"

for e in "${burstExe}" "${sasaExe}";
do
  if [ ! -x "$e" ]
  then
    echo "Error: setup error $e" 1>&2;
    exit 1;
  fi
done


usage(){
  echo "
================================================================
 Usage:  $0 -p <trajector pdbs> 
           [-c <chain id>] 
           [-l <file containing residue list>]
           [-s] 
================================================================
  l: residue ids of the Sasa of interested residues 
  s: Scwrl4 for rotamer adjustment
" 1>&2;
}

cleanup(){ rm -fr $hdir;}

halt(){ cleanup; usage; exit 1;}

hdir="/tmp/temporary_$$"

mkdir -p "$hdir"

while getopts ":p:c:l:s" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})"
        ;;
     c)
        chain=${OPTARG};
        ;;
     l)
        listfile="$(readlink -f ${OPTARG})";
        ;;
     s)
        rotamer=1;
        ;;
     *)
        halt;
        ;;
  esac
done
shift $((OPTIND-1))


if [ -z "${pdbfile}" ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ]
then
  halt;
fi

if [ ! -z $listfile ]
then
  if [ ! -f "${listfile}" ] || [ $(cat "$listfile" | wc -l) -ne $( grep -P '^\d+$' "$listfile" | wc -l) ]
  then
    halt;
  fi
fi

if [ $(grep -E "^MODEL" "${pdbfile}" | wc -l) -le 1 ]
then
  halt;
fi


"${burstExe}" -p "${pdbfile}" -o "$hdir" -e "snapshots" 

ls $hdir  | sort -t_ -k2n | while read pdb;
do
  modeId=$(basename $pdb | sed -e 's/snapshots//' -e 's/.pdb//' -e 's/_//g');
  pdbx="${hdir}/$pdb"
  cmd="${sasaExe} -p ${pdbx}"
  if [ ! -z $chain ]
  then
    cmd="$cmd -c $chain"
  fi
  
  if [ ! -z $listfile ]
  then
    cmd="$cmd -l $listfile"
  fi
  
  if [ ! -z $rotamer ]
  then
    cmd="$cmd -s"
  fi
  printf '%d,%.2f\n' $modeId $($cmd)
done

cleanup;
exit 0;
