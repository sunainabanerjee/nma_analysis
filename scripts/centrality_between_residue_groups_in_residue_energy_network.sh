#! /bin/bash

exedir=$(dirname $(readlink -f $0))
resCont="${exedir}/energy_based_weighted_contact.sh"

if [ ! -x "${resCont}" ]
then
  echo "Error: Setup issue" 1>&2;
  exit 1
fi

usage(){ 
 echo "
===============================================================================
Usage: $0 -p <pdb model file> 
         [-c <chain>] 
         [-r <contact radius>] 
         [-a <site1 resids file>] 
         [-b <site2 resids file>] 
         [-s] 
         [-m]
===============================================================================
c: one letter chain identifier
r: default contact radius 4.5
a: file containing list of residue ids at site - 1
b: file containing list of residue ids at site - 2
s: enable rotamer readjustment using scwrl
m: consider mainchain-mainchain interaction in the energy
" 1>&2;
}


rotamer=""
mainchain=""

while getopts ":p:c:a:b:r:sm" o;
do
  case "${o}" in
     p)
        pdbfile=$(readlink -f ${OPTARG})
        ;;
     c)
        chain=${OPTARG};
        ;;
     a)
        site1=$(readlink -f ${OPTARG});
        ;;
     b) 
        site2=$(readlink -f ${OPTARG});
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

if [ -z "${site1}" ] || [ ! -f "${site1}" ]
then
  usage;
  exit 1;
fi

if [ -z "${site2}" ] || [ ! -f "${site2}" ]
then
  usage;
  exit 1;
fi

if [ $(echo -n "$chain" | wc -c) -ne 1 ]
then
  usage;
  exit 1;
fi

if [ -z $radius ] || [ "$(echo $radius | grep -P '^\d+\.?\d*$')" != "$radius" ]
then
   radius=4.5
fi

ufile1="${site1}"
ufile2="${site2}"

target_file1="$(tempfile)"
target_file2="$(tempfile)"
contact_file="$(tempfile)"
tfile="$(tempfile)"
qfile="$(tempfile)"
nfile="$(tempfile)"

find_neighbor_pdb -f "${pdbfile}" -c ${chain} -e 2>/dev/null > $nfile

"${resCont}" -p ${nfile} -c ${chain} -r ${radius} ${rotamer} ${mainchain} > ${contact_file}

awk -vfile="${ufile1}" -vchain=${chain} 'BEGIN{ 
                                    while(getline line < file) 
                                       map[line] = 1;
                                 } /^ATOM/ && substr($0,22,1) == chain {
                                   resname= substr($0,18,3); 
                                   residue=substr($0,23,4); 
                                   gsub(" ","",residue); 
                                   if(map[residue] == 1 && found[residue] == 0){
                                        print resname""residue;
                                        found[residue] = 1;
                                   }
                                 }' "${pdbfile}" | uniq > ${target_file1}

awk -vfile="${ufile2}" -vchain=${chain} 'BEGIN{ 
                                    while(getline line < file) 
                                       map[line] = 1;
                                 } /^ATOM/ && substr($0,22,1) == chain {
                                   resname= substr($0,18,3); 
                                   residue=substr($0,23,4); 
                                   gsub(" ","",residue); 
                                   if(map[residue] == 1 && found[residue] == 0){
                                        print resname""residue;
                                        found[residue] = 1;
                                   }
                                 }' "${pdbfile}" | uniq > ${target_file2}


echo "
library(igraph)

x1 <- read.csv( \"${target_file1}\", header=F, sep=\",\" )
x2 <- read.csv( \"${target_file2}\", header=F, sep=\",\" )
y  <- read.csv( \"${contact_file}\", header=F, sep=\",\" )
el <- matrix(nrow = length(y\$V1), ncol = 2, )
el[,1] <- as.character(y\$V1)
el[,2] <- as.character(y\$V2)
g  <- graph_from_edgelist(as.matrix(el), directed = F)
mx <- max(y\$V3)
E(g)\$weight <- (mx + 1) - y\$V3

vertices <- names(V(g))

site1 <- c()
for( n in as.character(x1\$V1) ){
  if( n %in% vertices ){ site1 <- append(site1, n) }
}

site2 <- c()
for(n in as.character(x2\$V1)){
  if( n %in% vertices){ site2 <- append(site2,n) }
}

s  <- c()
for( n in site1 ){
  result <- shortest_paths( g, n, site2)\$vpath
  for(p in result){
    k <- unlist(p);
    if(length(k)>2){
      s <- append(s, names(k[2:(length(k)-1)]))
    }
  }
}
tbl <- table(s)

for( n in names(tbl)){
  cat( paste(paste(n, tbl[n], sep=\",\"),\"\\n\",sep=\"\") )
}

" > ${tfile}

Rscript "${tfile}" > "${qfile}"

sort -t, -k2,2nr $qfile

rm -f ${nfile} ${target_file1} ${taget_file2} ${contact_file} ${tfile} ${qfile} ${nfile}

