#! /bin/bash

exedir=$(dirname $(readlink -f $0));

modelBurstExe="${exedir}/burst_model.sh"

if [ ! -x "${modelBurstExe}" ]
then
  echo "Error: setup error [missing ${modelBurstExe}]" 1>&2;
  exit 1;
fi

usage(){ echo "
===================================================================================
Usage: $0 -p <pdb file> 
         [-c <chain id>] 
         [-o <out directory>] 
         [-e <extension >] 
	 [-k <number of modes to keep>]
	 [-b]
===================================================================================
default extension is 'model'
b: is tag for trigger create each mode model snapshot to burst into pdb instances.
" 1>&2; }

cleanup(){ rm -f $tfile $rfile; }

halt(){ cleanup; usage; exit 1; }

extension='model'

while getopts ":p:c:o:e:k:b" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
        ;;
     c)
        chain=${OPTARG};
        ;;
     o)
        outdir="${OPTARG}";
        ;;
     e)
        extension="${OPTARG}";
        ;;
     b)
        perform_burst=1;
        ;;
     k)
        keep=${OPTARG};
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1));


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

if [ -z $outdir ]
then
  outdir=$(basename "${pdbfile}" | sed -e 's/\.pdb$//');
fi

if [ -z $keep ] || [ "$(echo $keep|grep -P '^\d+$')" == "$keep" ]
then 
  keep=12
fi

if [ ! -d "${outdir}" ]
then
  mkdir -p "${outdir}"
fi

tfile=$(tempfile);
rfile=$(tempfile);

cd "${outdir}";

find_neighbor_pdb -f "${pdbfile}" -c ${chain} -e > $tfile

echo "
library(bio3d)
pdb <- read.pdb(\"${tfile}\")
modes <- aanma(pdb,outmodes=\"noh\", keep=$keep)
for(i in seq(7, length(modes))){
  file.name <- sprintf(\"mode_%d.pdb\", i)
  mktrj.nma(modes, mode=i, file=file.name, pdb=pdb)
}
" > "${rfile}"

Rscript "${rfile}"

if [ ! -z $perform_burst ]
then
   for pdbName in mode_*.pdb;
   do
     outDir=$(basename "${pdbName}")
     outDir=${outDir%.pdb}

     if [ ! -d "${outDir}" ]
     then
       mkdir -p "${outDir}"
     fi
  
     "${modelBurstExe}" -p "${pdbName}" -o "${outDir}" -e "${extension}"
   done
fi

cleanup;
