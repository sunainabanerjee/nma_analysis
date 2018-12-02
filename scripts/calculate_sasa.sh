#! /bin/bash

freeSasa=$(readlink -f $(which freesasa))
scwrl=$(readlink -f $(which Scwrl4))
findnbr=$(readlink -f $(which find_neighbor_pdb))

for e in "${freeSasa}" "${scwrl}" "${findnbr}";
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
 Usage:  $0 -p <pdb file> 
           [-c <chain id>] 
           [-l <file containing residue list>]
           [-s] 
================================================================
  l: residue ids of the Sasa of interested residues 
  s: Scwrl4 for rotamer adjustment
" 1>&2;
}

cleanup(){ rm -f $tfile $efile $qfile $zfile;}

halt(){ cleanup; usage; exit 1;}


tfile=$(tempfile)
efile=$(tempfile)
qfile=$(tempfile)
zfile=$(tempfile)


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


if [ -z "${chain}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain=$(find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi

if [ $(echo -n "$chain" | wc -c) -ne 1 ]
then
  halt;
fi

if [ $(grep -E "^MODEL" "$pdbfile" | wc -l) -gt 1 ]
then
  halt;
fi

if [ ! -z $listfile ]
then
  if [ ! -f "${listfile}" ] || [ $(cat "$listfile" | wc -l) -ne $( grep -P '^\d+$' "$listfile" | wc -l) ]
  then
    halt;
  fi
else
  awk -vchain=$chain '/^ATOM/ && substr($0,22,1) == chain { 
                          resid=sprintf("%d", substr($0,23,4)); 
                          print resid;}' "${pdbfile}"  | sort -nu > $zfile
  listfile=$zfile
fi


"${findnbr}" -f "${pdbfile}" -c $chain -e > $tfile 2>/dev/null

if [ -z $rotamer ]
then
  cp "$tfile" "$efile"
else
  "${scwrl}" -i "${tfile}" -o "${efile}" >/dev/null
fi

"${freeSasa}" "${efile}"  --format=rsa | awk '/^RES/{ 
                                                  printf("%s,%d,%.2f\n",$2,$4,$5);
                                             }' | awk -vfile="$listfile" -F, 'BEGIN{ 
                                                                            while( getline line < file ){
                                                                              residues[line]=1;
                                                                            }
                                                                           }NF==3 && residues[$2] == 1{
                                                                              s += $3;
                                                                           }END{
                                                                             print s;
                                                                           }'

cleanup;
exit 0;
