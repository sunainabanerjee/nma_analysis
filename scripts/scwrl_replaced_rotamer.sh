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
================================================================================
" 1>&2; }

cleanup(){ rm -fr $zfile $outdir; }

halt(){ usage; cleanup; exit 1; }


extension="model";

while getopts ":p:" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
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

zfile=$(tempfile)
outdir="/tmp/tmpdir_$$"

if [ ! -d $outdir ]
then
   mkdir -p "$outdir" || halt;
fi

grep '^MODEL' "${pdbfile}"  | while read line; 
do
   snapShot=$(echo $line | awk '{print $2}');
   ofile="${outdir}/${extension}_${snapShot}.pdb"
   sed -n "/$line/,/ENDMDL/p"  "${pdbfile}" > "${ofile}";
   "${scwrlExe}" -i "${ofile}" -o "${zfile}" >/dev/null
   echo "$line"   > "${ofile}"
   cat "${zfile}" >> "${ofile}"
   echo "ENDMDL"  >> "${ofile}"
done

ls ${outdir}/${extension}_*.pdb | sort -t"_" -k3n | xargs -i cat {}

cleanup;
exit 0;
