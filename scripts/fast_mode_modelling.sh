#! /bin/bash

exeDir=$(dirname $(readlink -f $0))
scwrlExe=$(readlink -f $(which Scwrl4))
pulchraExe=$(readlink -f $(which pulchra))
burstExe="${exeDir}/burst_model.sh"

for e in "${scwrlExe}" "${pulchraExe}" "${burstExe}"
do
  if [ ! -x "${e}" ]
  then
    echo "Error: executable not found $e " 1>&2;
    exit 1;
  fi
done


usage(){ echo "
================================================================================
Usage: $0 -p <pdb model file> 
          -o <model out directory> 
          [-x <outmodel extension>]
          [-s <start mode>]
          [-e <end mode>]
          [-f <force field>]
          [-m]
          [-S]
================================================================================
p :  expects pdb with single pdb Chain only
e :  default extension is \"model\", in output directory 
     all model pdbs will be named as model_<model id>.pdb
x :  default extension is 'mode'
s :  default start mode is 7
e :  default end mode is 12
f :  forcefield to be used, default value is 'sdenm', valid 
     values ('calpha', 'anm', 'reach', 'sdenm' )
m :  perform basic backbone modelling (default is off)
S :  switch on scwrl rotamer readjustment
" 1>&2; }

cleanup(){ rm -fr $tdir $tfile $rfile; }

halt(){ usage; cleanup; exit 1; }

extension="mode";
smode=7
emode=12
ff=sdenm

while getopts ":p:e:o:s:e:f:x:Sm" o;
do
  case "${o}" in
     p)
        inputPdb="$(readlink -f ${OPTARG})";
        ;;
     o)
        outdir="${OPTARG}";
        ;;
     x)
        extension="${OPTARG}";
        ;;
     s)
        smode=${OPTARG};
        ;;
     e) 
        emode=${OPTARG};
        ;;
     f)
        ff=${OPTARG};
        ;;
     S)
        rotamer=1;
        ;;
     m)
        perform=1;
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1));

if [ -z ${inputPdb} ] || [ ! -f "${inputPdb}" ] || [ $(awk '/^ATOM/' "${inputPdb}"| wc -l) -lt 4 ]
then
   halt;
fi

if [ -z ${outdir} ] || [ ! -d "${outdir}" ]
then
   halt;
fi

tfile=$(tempfile)
rfile=$(tempfile)
tdir=/tmp/tempdir_$$

if [ -d $tdir ]
then
  rm -fr $tdir
fi
mkdir -p $tdir

if [ ! -d "${outdir}" ]
then
  mkdir "${outdir}"
fi


echo "
 library(bio3d)
 pdb <- read.pdb('${inputPdb}')
 modes <- nma(pdb, ff='$ff')
 for( i in seq($smode,$emode) ){
   fname <- sprintf(\"%s/%s_catraj_%d.pdb\", \"${outdir}\", \"${extension}\", i)
   mktrj(modes, mode=i, pdb=pdb, rock=TRUE, step=.5, file=fname)
 }
" > $rfile

Rscript $rfile


if [ ! -z $perform ]
then
 for f in $outdir/${extension}_catraj_*.pdb
 do
   modeId=$(basename $f | sed -e 's/mode_catraj_\([0-9]\+\).pdb/\1/')
   burstDir=$tdir/burst_$$
   
   if [ -d $burstDir ]
   then
     rm -fr $burstDir
   fi
   mkdir -p $burstDir

   $burstExe -p "$f" -o "$burstDir" -e "model"

   ls $burstDir/model_*.pdb | grep -P 'model_\d+.pdb' | 
       sed -e 's/\(.*model_\([0-9]\+\).pdb\)$/\2 \1/' | 
                                            sort -k1n | 
                                      awk '{print $2}'| while read j;
   do
     modelId=$(basename $j | sed -e 's/model_\([0-9]\+\).pdb/\1/')
     awk '/^ATOM/ {
        x=sprintf("%.3f",substr($0,30,8)); 
        y=sprintf("%.3f",substr($0,39,8)); 
        z=sprintf("%.3f", substr($0,47,8)); 
        print x,y,z;
     }' $j   > $tfile
   
     awk -vfile=$tfile 'BEGIN{ 
        i=1; 
        while( getline line < file){ 
        split(line,arr," "); 
        if( length(arr)==3 ){ 
           X[i]=arr[1]; 
           Y[i]=arr[2]; 
           Z[i]=arr[3];
           i++;
        }
      } 
     }
     /^ATOM/ && substr($0,14,3) ~ /^CA/ && substr($0,17,1) ~ /(A| )/{ 
       counter++;
       printf("%s%8.3f%8.3f%8.3f%s\n",substr($0,1,30),X[counter],Y[counter],Z[counter],substr($0,55)); 
     }' "${inputPdb}" > $j
     pulchra $j 2>&1 1>/dev/null
     if [ ! -z $rotamer ]
     then
        Scwrl4 -i ${j%.pdb}.rebuilt.pdb -o $j 2>&1 1>/dev/null
     else
        cp ${j%.pdb}.rebuilt.pdb $j
     fi
     awk -vmodel=$modelId 'BEGIN{ printf("MODEL     %4d\n",model);} /^ATOM/{ print } END{ printf("ENDMDL\n")}' $j
   done > $outdir/${extension}_${modeId}.pdb
 done
fi

rm -fr $tdir $tfile $rfile
