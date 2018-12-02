#! /bin/bash


usage(){ echo "Usage: $0 -p <pdb file> [-c <chain id>] -i <residue mapfile>" 1>&2; }



while getopts ":p:c:i:" o;
do
  case "${o}" in 
     p)
        pdbfile="$(readlink -f ${OPTARG})"
        ;;
     c)
        chain=${OPTARG};
        ;;
     i)
        residuemap="$(readlink -f ${OPTARG})";
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${pdbfile}" ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ]
then
  usage;
  exit 1;
fi

if [ -z "${residuemap}" ]  || [ ! -f "${residuemap}" ]
then
  usage;
  exit 1;
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
  usage;
  exit 1;
fi

nfile=$(tempfile)

find_neighbor_pdb -f "${pdbfile}" -c "${chain}" -e 2>/dev/null >$nfile

awk -vsfile="${residuemap}" 'BEGIN{
                          while(getline line < sfile){
                            split(line, arr, ",");
                            mx = 0;
                            if( length(arr) >= 2 ){
                              resId=substr(arr[1],4);
                              entropy[resId]=arr[2];
                            }
                          }
                       }
                       /^ATOM/{
                          residue=sprintf("%d", substr($0,23,4));
                          color = entropy[residue];
                          printf "%s%6.2f%s\n",substr($0,1,60),color,substr($0,67);
                       }' "${nfile}"

rm -f $nfile

