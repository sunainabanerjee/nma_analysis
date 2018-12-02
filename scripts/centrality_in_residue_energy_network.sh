#! /bin/bash

exedir=$(dirname $(readlink -f $0))
resCont="${exedir}/energy_based_weighted_contact.sh"
unwtCentExe="${exedir}/weighted_centrality.sh"


for e in "${resCont}" "${unwtCentExe}";
do
  if [ ! -x "${e}" ]
  then
    echo "Error: Setup issue [Missing executable $e]" 1>&2;
    exit 1
  fi
done

usage(){ 
 echo "
------------------------------------------------------------------
Usage: $0 -p <pdb file> [-c <chain>] [ -k <centrality type>] 
          [-r <contact radius>] [-s] [-m]
------------------------------------------------------------------
k: centrality measure 'betweenness', 'pagerank', 'closeness', 'eigen',
   'subgraph'
   default 'betweenness'
r: default radius cutoff is 4.5
s: tag switches the rotamer correction (using scrwl4)
m: switched on main chain main chain interaction
" 1>&2;
}

cleanup(){ rm -f ${contact_file} ${nfile}; }
rotamer=""
mainchain=""

while getopts ":p:c:k:r:sm" o;
do
  case "${o}" in
     p)
        pdbfile=$(readlink -f ${OPTARG})
        ;;
     c)
        chain=${OPTARG};
        ;;
     k)
        centrality=${OPTARG};
        ;;
     r)
        radius=${OPTARG};
        ;;
     s)
        rotamer="-s";
        ;;
     m)
        mainchain="-m";
        ;;
     *)
        usage;
        exit 1;
        ;;
  esac
done
shift $((OPTIND-1))

if [ -z $centrality ] || [ $(echo "$centrality" | grep -E '^(betweenness|eigen|pagerank|closeness|subgraph)$') != "${centrality}" ]
then
  centrality='betweenness'
fi

if [ -z "${pdbfile}" ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ]
then
  usage;
  exit 1;
fi

if [ -z $radius ] || [ "$(echo $radius | grep -P '^\d+\.\d+$')" != "${radius}" ]
then 
  radius=4.5
fi

if [ -z "${chain}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain=$(find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi


contact_file="$(tempfile)"
nfile="$(tempfile)"

"${resCont}" -p ${pdbfile} -c ${chain} -r ${radius} ${rotamer} ${mainchain} > ${contact_file}

"${unwtCentExe}" -f "${contact_file}" -k "${centrality}"

cleanup;
exit 0;

