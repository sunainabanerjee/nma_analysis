#! /bin/bash


usage(){ echo "Usage: $0 -p <pdb file> [-c <chain id>]" 1>&2; }



while getopts ":p:c:" o;
do
  case "${o}" in 
     p)
        pdbfile=$(readlink -f ${OPTARG})
        ;;
     c)
        chain=${OPTARG};
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


awk -vchain=$chain '/^ATOM/ && substr($0,22,1) == chain {
                       resId=sprintf("%d", substr($0,23,4));
                       resName=substr($0,18,3); 
                       restag=sprintf( "%s%d", resName, resId);
                       if( seen[restag] == 0)
                          print restag;
                       seen[restag] = 1;
                   }' "${pdbfile}" | awk '{printf("%d,%s\n",NR,$1);}'
