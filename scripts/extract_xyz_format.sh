#! /bin/bash

usage(){ echo "Usage: $0 -p <pdb file> [-c <chain id>]" 1>&2; }

cleanup(){ rm -f $tfile; }

halt(){ cleanup; usage; exit 1; }


exedir="$(dirname $(readlink -f "$0"))"
tfile=$(tempfile)


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

find_neighbor_pdb -f "${pdbfile}" -c "$chain" -e > $tfile

awk '/^ATOM/ && substr($0,17,1) ~ /( |A)/{
   counter++;
   printf("%3s %4d %4s %4d %8.3f %8.3f %8.3f \n",
                   substr($0,18,3),
                   substr($0,23,4),
                   substr($0,13,4), 
                   counter,
                   substr($0,31,8), 
                   substr($0,39,8),
                   substr($0,47,8)); }' "$tfile"

cleanup;
