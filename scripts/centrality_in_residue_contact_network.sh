#! /bin/bash

exedir=$(dirname $(readlink -f $0))
resCont="${exedir}/within_chain_residue_residue_network.sh"
unwtCentExe="${exedir}/unweighted_centrality.sh"


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
------------------------------------------------------------------
k: centrality measure 'betweenness', 'pagerank', 'closeness', 'eigen',
   'subgraph'
   default 'betweenness'
" 1>&2;
}

cleanup(){ rm -f ${contact_file} ${nfile}; }


while getopts ":p:c:k:" o;
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

if [ -z "${chain}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain=$(find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi


contact_file="$(tempfile)"
nfile="$(tempfile)"

find_neighbor_pdb -f "${pdbfile}" -c ${chain} -e 2>/dev/null > $nfile

"${resCont}" -p ${nfile} -c ${chain} > ${contact_file}

"${unwtCentExe}" -f "${contact_file}" -k "${centrality}"



