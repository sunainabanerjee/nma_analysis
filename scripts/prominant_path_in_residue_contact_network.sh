#! /bin/bash

exedir=$(dirname $(readlink -f $0))
resCont=${exedir}/within_chain_residue_residue_network.sh

if [ ! -x "${resCont}" ]
then
  echo "Error: Setup issue" 1>&2;
  exit 1
fi

usage(){ echo "
------------------------------------------------------------------
Usage: $0 -p <pdb file> [-c <chain>] [-m <minimum path length>] 
------------------------------------------------------------------
default minimum path length 5
" 1>&2; }

cleanup(){ rm -f $contact_file $tfile $qfile $nfile $pfile;}

while getopts ":p:c:m:" o;
do
  case "${o}" in
     p)
        pdbfile=$(readlink -f ${OPTARG})
        ;;
     c)
        chain=${OPTARG};
        ;;
     m)
        min_pathlen=${OPTARG};
        ;;
     *)
        usage;
        exit 1;
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

if [ -z $min_pathlen ] || [ "$(echo "$min_pathlen" | grep -P '^\d+$')" != "${min_pathlen}" ]
then
  min_pathlen=5
fi

contact_file="$(tempfile)"
tfile="$(tempfile)"
qfile="$(tempfile)"
nfile="$(tempfile)"
pfile="$(tempfile)"

find_neighbor_pdb -f "${pdbfile}" -c ${chain} -e 2>/dev/null > $nfile

"${resCont}" -p ${nfile} -c ${chain} > ${contact_file}

echo "
library(igraph)

y  <- read.csv( \"${contact_file}\", header=F, sep=\",\" )
g  <- graph_from_edgelist(as.matrix(y), directed = F)
x  <- components(g)

m  <-  ${min_pathlen}

for( i in seq(length(x\$csize) ) ){
  gn <- names( x\$membership[x\$membership == i] )
  n <- length(gn)
  for( i1 in seq(n-1) ){
    for( i2 in seq(i1+1, n)){
       p <- names( shortest_paths(g, gn[i1], gn[i2])\$vpath[[1]] )
       if( length(p) > m + 2 ){
          for(i in seq(length(p)-1)){ 
               a <- as.numeric(substring(p[i],  4)); 
               b <- as.numeric(substring(p[i+1],4)); 
               tag <- ifelse( a > b , sprintf(\"%s:%s\", p[i], p[i+1]), 
                              sprintf(\"%s:%s\",p[i+1], p[i]) ); 
               cat(paste(tag,'\\n'))
          }
       }
    }
  }
}
" > $tfile

Rscript "${tfile}" > "${qfile}"

sort $qfile | uniq -c | sort -k1,1nr | awk 'NF==2' > $tfile

sort -k1,1nr $tfile | awk  'NF==2 { 
                               if(NR==1) 
                                 max = $1 + 1; 
                               split($2,arr, ":"); 
                               printf("%s,%s,%d\n", arr[1],arr[2],max - $1);
                            }' > $qfile

max=$(tail -1 $qfile | cut -d, -f3)

echo "
library(igraph)

y  <- read.csv( \"$qfile\", header=F, sep=\",\" )
el <- matrix(nrow = length(y\$V1), ncol = 2, )
el[,1] <- as.character(y\$V1)
el[,2] <- as.character(y\$V2)
g  <- graph_from_edgelist(as.matrix(el), directed = F)
E(g)\$weight <- y\$V3
g.tree <- mst(g)
for( e in E(g.tree)){ 
  verts <- get.edge(g.tree,e); 
  score <- $max - E(g.tree)\$weight[e]
  cat( paste( paste( names(V(g.tree))[verts[1]], names(V(g.tree))[verts[2]], score, sep=\",\" ), \"\\n\", sep=\"\") );
}
" > $tfile

Rscript "${tfile}"

cleanup;

