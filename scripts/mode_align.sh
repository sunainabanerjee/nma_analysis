#! /bin/bash

error_msg() {
  echo "Error: $*" 1>&2;
}

symmetric_diff () {
  tfile1=$(tempfile);
  tfile2=$(tempfile);
  awk -F, -vOFS=":" '{ print $1,$2 }' "$1" > $tfile1
  awk -F, -vOFS=":" '{ print $1,$2 }' "$2" > $tfile2
  cat $tfile1 $tfile2 | sort | uniq -c | awk '$1 == 1' | wc -l
  rm -f $tfile1 $tfile2
}


findNbr=$(readlink -f $(which find_neighbor_pdb))
exedir=$(dirname $(readlink -f $0))
modeAna=${exedir}/fast_mode_modelling.sh
allPairExe=${exedir}/allpair_mobility

for exe in "${findNbr}" "${modeAna}" "${allPairExe}"
do
  if [ -f "${exe}" ] && [ ! -x "${exe}" ]
  then
    error_msg "Missing executable $exe" ;
    exit 1;
  fi
done


usage(){ echo "
================================================================================
Usage: $0 -p <pdb file>
          -o < output directory > 
          [ -e < end mode > ]
          [ -f < force field > ]
================================================================================
e:  end mode (default value: 25)
f:  force field to be used for the NMA analysis [calpha, sdenm, reach] (default: sdenm)
o:  output directory ( default: <input>_mutate.pdb )
" 1>&2; }

cleanup(){ rm -f $tfile $rfile $qfile; }

halt(){ usage; cleanup; exit 1; }

smode=7

while getopts ":p:e:o:f:" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
        ;;
     o)
        outdir="$(readlink -f ${OPTARG})";
        ;;
     e)
        emode=${OPTARG};
        ;;
     f)
        ff=${OPTARG};
        ;;
     *)
        usage;
        exit 0;
        ;;
 esac
done
shift $((OPTIND-1));


if [ -z ${pdbfile} ] || [ ! -f "${pdbfile}" ] ||  [ "$(echo $pdbfile | grep -E ".pdb$")" != "$pdbfile" ]
then
   halt;
fi

if [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ] 
then 
  error_msg "Improper pdb file content [$pdbfile]";
  exit 1;
fi

if [ $( $findNbr -f "$pdbfile" -C | wc -l) -ne 1 ] 
then
  error_msg "Multiple chains in the pdb file [$pdbfile]" ;
  exit 1;
fi

if [ -z ${outdir} ] 
then
   halt;
fi

if [ -z $ff ] || [ "$( echo $ff | grep -E "^(sdenm|calpha|reach)")" != "$ff" ]
then
  ff=sdenm
fi

if [ -z $emode ] || [ "$(echo $emode | grep -E "^\d+$" )" != "$emode" ]
then
  emode=25
fi

if [ ! -d $outdir ]
then
  mkdir -p $outdir
fi

tfile=$(tempfile)
rfile=$(tempfile)
qfile=$(tempfile)

x=$(find_neighbor_pdb -f "$pdbfile" -C | awk -F: '{print $2}')
x=$((x*(x-1)/2))
x=$((x/75))

$modeAna -p "${pdbfile}" -o "${outdir}" -f $ff -s $smode -e $emode

for f in $outdir/mode_catraj_*.pdb
do
  if [ -f "$f" ] && [ ! -f "${f%pdb}scr" ]
  then
   awk '/^MODEL/{ 
        modelId=$2; 
     } 
     /^ATOM/ && substr($0,14,3) ~ /CA / { 
        resId=sprintf("%d", substr($0,23,4)); 
        name=substr($0,18,3);
        x = sprintf("%.3f", substr($0,30,8)); 
        y = sprintf("%.3f", substr($0,39,8)); 
        z = sprintf("%.3f", substr($0,47,8)); 
        X[resId,modelId] = x; 
        Y[resId,modelId] = y; 
        Z[resId,modelId] = z; 
        residues[resId] = 1; 
        model[modelId] = 1;
        resName[resId] = name;
     } END{
          for( m in model){
            for( r in  residues){
               printf("%d,%d,%.3f,%.3f,%.3f\n", m, r, X[r,m], Y[r,m], Z[r,m] );
            }
          } 
     }' "$f" > $qfile
     $allPairExe  "$qfile"  | sort -t, -k5 -nr | head -${x} > ${f%pdb}scr
  fi
done

exit 1;

for i in ${outdir}/mode_catraj*.scr; 
do 
   iId=$(basename $i | sed -e 's/mode_catraj_\([0-9]\+\).scr/\1/')
   for j in ${outdir}/mode_catraj*.scr; 
   do 
      jId=$(basename $j | sed -e 's/mode_catraj_\([0-9]\+\).scr/\1/')
      if [[ "$iId" -lt "$jId" ]]; 
      then 
        echo $iId,$jId,$(symmetric_diff "$i" "$j"); 
      fi 
   done 
done > $tfile


awk -F, -vnorm=$x -vstart=$smode -vend=$emode 'NF==3{ 
                                score[$1,$2] = sprintf( "%f", $3 * 1.0/norm); 
                                score[$2, $1] = score[$1, $2];
                             }END{
                               for( i = start; i <= end; ++i ){
                                 for( j = start; j <= end; ++j ){
                                    if( j != start )
                                       printf(",");
                                    printf("%.3f", score[i,j]);
                                 }
                                 printf("\n");
                               }
                             }' $tfile > "${outdir}/mode_pair_score.csv"

echo "
library(ggplot2)
library(ggdendro)

x <- read.csv( \"$tfile\", head=F )
x.names <- sort(unique(c(x\$V1, x\$V2)))
x.dist <- matrix(0, length(x.names), length(x.names))
dimnames(x.dist) <- list(x.names, x.names)
x.ind <- rbind(cbind(match(x\$V1, x.names), 
                     match(x\$V2, x.names)),
                     +      cbind(match(x\$V2, x.names), match(x\$V1, x.names)))
x.dist[x.ind] <- rep(x\$V3, 2)

hc <- hc <- hclust(as.dist(x.dist))
pdf( \"${outdir}/mode_motion_comparison.pdf\")
ggdendrogram(hc, size=2)
dev.off()

x.dist <- x.dist / (2. * $x)
write.table( x.dist , file=\"${outdir}/mode_motion_comparison.csv\", sep=',' , col.names=TRUE, row.names=TRUE)
" > $rfile

Rscript $rfile

cleanup;
