#! /bin/bash

scwrlExe=$(readlink -f $(which Scwrl4))

for e in "${scwrlExe}"
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
          [-e <outmodel extension>]
          [-s]
================================================================================
e : default extension is \"model\", in output directory 
    all model pdbs will be named as model_<model id>.pdb
s:  switch on scwrl rotamer readjustment
" 1>&2; }

cleanup(){ rm -f $zfile; }

halt(){ usage; cleanup; exit 1; }


extension="model";

while getopts ":p:e:o:s" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
        ;;
     o)
        outdir="${OPTARG}";
        ;;
     e)
        extension="${OPTARG}";
        ;;
     s)
        rotamer=1;
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1));


if [ -z ${pdbfile} ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ]
then
   halt;
fi

if [ -z ${outdir} ] || [ ! -d "${outdir}" ]
then
   halt;
fi

zfile=$(tempfile)

grep '^MODEL' "${pdbfile}"  | while read line; 
do
   snapShot=$(echo $line | awk '{print $2}');
   ofile="${outdir}/${extension}_${snapShot}.pdb"
   sed -n "/$line/,/ENDMDL/p"  "${pdbfile}" > "${ofile}";
   if [ ! -z $rotamer ]
   then
     "${scwrlExe}" -i "${ofile}" -o "${zfile}" >/dev/null
     mv "${zfile}" "${ofile}"
   fi
done

cleanup;
exit 0;
