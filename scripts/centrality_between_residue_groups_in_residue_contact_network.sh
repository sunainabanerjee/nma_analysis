#! /bin/bash

exedir=$(dirname $(readlink -f $0))
resCont=${exedir}/within_chain_residue_residue_network.sh

if [ ! -x "${resCont}" ]
then
  echo "Error: Setup issue" 1>&2;
  exit 1
fi

usage(){ 
 echo "
==============================================================
Usage: $0 
          -p <pdb file> 
         [-c <chain>] 
         [-a <site1>] 
         [-b <site2>]
==============================================================
c: one letter chain identifier of the molecule
   in question. This is optional for pdb file with 
   only single chain.
a: site residue ids at site 1
b: site residue ids at site 2
" 1>&2;
}


while getopts ":p:c:a:b:" o;
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


ufile1="${site1}"
ufile2="${site2}"

target_file1="$(tempfile)"
target_file2="$(tempfile)"
contact_file="$(tempfile)"
tfile="$(tempfile)"
qfile="$(tempfile)"
nfile="$(tempfile)"

find_neighbor_pdb -f "${pdbfile}" -c ${chain} -e 2>/dev/null > $nfile

"${resCont}" -p ${nfile} -c ${chain} > ${contact_file}

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
g  <- graph_from_edgelist(as.matrix(y), directed = F)

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

cat $qfile

rm -f ${nfile} ${target_file1} ${taget_file2} ${contact_file} ${tfile} ${qfile} ${nfile}

