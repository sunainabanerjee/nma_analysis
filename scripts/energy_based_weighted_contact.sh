#! /bin/bash

exedir=$(dirname $(readlink -f $0));

modelBurstExe="${exedir}/burst_model.sh"
nbrEnergyExe="${exedir}/find_interacting_residue_energy"
xyzExtractor="${exedir}/extract_xyz_format.sh"
scwrlExe="$(which Scwrl4)"


for e in "${modelBurstExe}" "${nbrEnergyExe}" "${xyzExtractor}" "${scwrlExe}";
do
  if [ ! -x "${e}" ] 
  then
    echo "Error: setup error [missing ${e}]" 1>&2;
    exit 1;
  fi
done

usage(){ echo "
============================================================================
Usage: $0 -p <pdb model file> [-c <chain>] [-r <interaction influence radius>] [-s] [-m]
============================================================================
r: default interaction influence radius is 4.5 A
s: switch on scwrl side chain replacement or not
m: include main chain mainchain interaction into consideration
" 1>&2; }

cleanup(){ rm -fr $tdir $tfile $rfile $zfile; }

halt(){ cleanup; usage; exit 1; }

extension='model'
mainchain=""

while getopts ":p:r:c:sm" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
        ;;
     r)
        radius=${OPTARG};
        ;;
     s)
        rotamer=1;
        ;;
     m)
        mainchain="-m";
        ;;
     c)
        chain=${OPTARG};
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1));


if [ -z "${pdbfile}" ] || [ ! -f "${pdbfile}" ] || [ $(grep '^ATOM  ' "${pdbfile}" | wc -l) -lt 4 ]
then
   halt;
fi

if [ $(grep '^MODEL ' "${pdbfile}" | wc -l) -gt 1 ]
then
  model=1
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


if [ -z ${radius} ] || [ "$(echo $radius | grep -P '^\d+\.?\d*')" != "${radius}" ]
then
   radius=4.5;
fi

tdir="/tmp/tdir_$$"
tfile="$(tempfile)"
rfile="$(tempfile)"
zfile="/tmp/temp_$$.pdb"

if [ ! -d "${tdir}" ]
then
  mkdir -p "${tdir}"
fi

if [ ! -z $model ]
then
  "${modelBurstExe}" -p "${pdbfile}" -o "${tdir}" -e "${extension}"
else
  find_neighbor_pdb -f "${pdbfile}" -c ${chain} -e > "${tdir}/${extension}_1.pdb" 
fi

nCounter=$(ls ${tdir}/${extension}_*.pdb | wc -l)

for snapShot in $tdir/${extension}_*.pdb;
do
  if [ ! -z $rotamer ]
  then
    "${scwrlExe}" -i "${snapShot}" -o "${zfile}" >/dev/null
    mv "${zfile}" "${snapShot}"
  fi
  "${xyzExtractor}" -p "${snapShot}" -c ${chain} > $tfile;
  "${nbrEnergyExe}" -f "$tfile" -r ${radius} ${mainchain} >> "${rfile}"
done

awk -F, -vn=$nCounter  'NF==3 {
                            residue[$1]=1; 
                            residue[$2]=1; 
                            interaction[$1,$2] += $3;
                            hit[$1,$2] += 1;
                        }END{ 
                          for( a in residue )
                           for( b in residue )
                             if( interaction[a,b] != 0 )
                               printf( "%s,%s,%.3f,%d\n", a, b, interaction[a,b]/n, hit[a,b] ); 
                        }' "${rfile}" 

cleanup;
exit 0;
