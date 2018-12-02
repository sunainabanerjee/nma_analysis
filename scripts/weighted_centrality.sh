#! /bin/bash

exedir=$(dirname $(readlink -f $0))

usage(){ 
 echo "
------------------------------------------------------------------
Usage: $0 -f <contact csv> [ -k <centrality type>]
------------------------------------------------------------------
k: centrality measure 'betweenness', 'pagerank', 'closeness', 'eigen',
   'subgraph'
   default 'betweenness'
" 1>&2;
}

cleanup(){ rm -f $tfile $qfile; }

halt(){ cleanup; usage; exit 1;}

check_network(){
  network="$1";
  n=$(awk -F, 'NF >= 3 && $1 ~ /^[A-Z]{3}\-?[0-9]+$/ && $2 ~ /^[A-Z]{3}\-?[0-9]+$/ && $3 ~ /^[0-9]*.?[0-9]+$/' "$network" | wc -l);
  N=$(cat "$network" | wc -l);
  if [ $n -eq $N ]
  then
     echo "TRUE";
  fi
}


while getopts ":f:k:" o;
do
  case "${o}" in
     f)
        contact_file=$(readlink -f ${OPTARG})
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

if [ -z $contact_file ] || [ ! -f "${contact_file}" ] || [ "$(check_network "${contact_file}")" != "TRUE" ]
then
  halt;
fi


if [ -z $centrality ] || [ $(echo "$centrality" | grep -E '^(betweenness|eigen|pagerank|closeness|subgraph)$') != "${centrality}" ]
then
  centrality='betweenness'
fi


tfile=$(tempfile)
qfile=$(tempfile)

echo "
library(igraph)

y  <- read.csv( \"${contact_file}\", header=F, sep=\",\" )
el <- matrix(nrow = length(y\$V1), ncol = 2, )
el[,1] <- as.character(y\$V1)
el[,2] <- as.character(y\$V2)
g  <- graph_from_edgelist(as.matrix(el), directed = F)
mx <- max(y\$V3)
E(g)\$weight <- (mx + 1) - y\$V3
" > ${tfile}

if [ "$centrality" == 'betweenness' ]
then
echo "
tbl <- betweenness(g, normalized=T) * 100
" >> ${tfile}
elif [ "$centrality" == 'closeness' ]
then
echo "
tbl <- closeness(g, normalized=T) * 100
" >> ${tfile}
elif [ "$centrality" == 'pagerank' ]
then
echo "
d <- page_rank(g)
tbl <- (d\$vector / max(d\$vector)) * 100
" >> ${tfile}
elif [ "$centrality" == 'subgraph' ]
then
echo "
d <- subgraph_centrality(g)
tbl <- (d / max(d)) * 100
" >> ${tfile}
else
echo "
tbl <- eigen_centrality(g, scale=T)\$vector * 100
" >> ${tfile}
fi

echo "
for( n in names(tbl)){
  cat( paste(paste(n, tbl[n], sep=\",\"),\"\\n\",sep=\"\") )
}
" >> ${tfile}

Rscript "${tfile}" > "${qfile}"

awk -F, '{printf("%s,%.3f\n",$1,$2);}' $qfile | sort -t, -k2,2nr

cleanup;
exit 0;
