#! /bin/bash


usage(){ echo "
================================================================================
Usage: $0 -f <contact csv file> [-c <minimum cluster size>]
================================================================================
c: default cluster cutoff size is 10
" 1>&2; }

cleanup(){ rm -f $tfile; }

halt(){ usage; cleanup; exit 1; }


while getopts ":f:c:" o;
do
  case "${o}" in
     f)
        energyContact="$(readlink -f ${OPTARG})";
        ;;
     c)
        contact="${OPTARG}";
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1));


if [ -z ${energyContact} ] || [ ! -f "${energyContact}" ] || [ $(awk -F, 'NF >= 3 && $1 ~ /^[A-Z]{3}-?[0-9]+$/ && $2 ~ /^[A-Z]{3}-?[0-9]+$/' "${energyContact}"| wc -l) -ne $(cat "${energyContact}" | wc -l) ]
then
   halt;
fi

if [ -z ${contact} ] || [ "$(echo $contact | grep -P "^\d+$")" != "${contact}" ]
then
  contact=10
fi


tfile=$(tempfile)

echo "
library(igraph)
fname <- \"${energyContact}\" 
size.cutoff <- $contact
y <- read.csv(fname, sep = \",\", header = F)
el <- matrix(nrow = length(y\$V1), ncol = 2 )
el[,1] <- as.character(y\$V1)
el[,2] <- as.character(y\$V2)
g  <- graph_from_edgelist(as.matrix(el), directed = F)
mx <- max(y\$V3)
E(g)\$weight <- (mx + 1) - y\$V3
x <- edge_betweenness(g, directed = F)

g2 <- graph_from_edgelist(as.matrix(el), directed=F)
mx <- max(x)
E(g2)\$weight <- (mx + 1) - x
g2.mst <-  mst(g2)
d <- density(y\$V3)
cutoff <- max(d\$x[which.max(d\$y)])
vlist1 <- c()
vlist2 <- c()
for( e in E(g2.mst)){ 
  verts <- ends(g2.mst,e)
  v1 <- verts[1]
  v2 <- verts[2]
  e  <- y\$V3[ (y\$V1 %in% c(v1,v2)) & (y\$V2 %in% c(v1,v2))]
  if( e > cutoff ){
    vlist1 <- append(vlist1, v1)
    vlist2 <- append(vlist2, v2)
  }
}
mtree <- graph_from_edgelist(as.matrix(data.frame(v1 = vlist1, v2=vlist2 )), directed = F)
mtree.comp <- components(mtree) 
component.ids <- which(mtree.comp\$csize >= size.cutoff)

for( u in component.ids){
  vnames <- names(mtree.comp\$membership[mtree.comp\$membership == u])
  g2.connected <- induced_subgraph(mtree, vnames)
  for( e in E(g2.connected)){
    verts <- ends(g2.connected, e)
    v1 <- verts[1]
    v2 <- verts[2]
    e  <- y\$V3[ (y\$V1 %in% c(v1,v2)) & (y\$V2 %in% c(v1,v2))]
    cat( paste(paste( v1, v2, e, u, sep=\",\" ),\"\n\", sep=\"\") )
  }
}
" > "${tfile}"

Rscript $tfile

cleanup;
exit 0;
