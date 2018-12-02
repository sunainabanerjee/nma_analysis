#! /bin/bash

findNbr=$(readlink -f $(which find_neighbor_pdb))
scwrl=$(readlink -f $(which Scwrl4))


for exe in "${findNbr}" "${scwrl}"
do
  if [ ! -x "${exe}" ]
  then
    echo "Error: missible executable ${exe}" 1>&2;
    exit 1;
  fi
done

usage(){ echo "
================================================================================
Usage: $0 -p <pdb file>
          [ -r < residue ID >  ]
          [ -s < sequence ID > ]
          [ -a < mutation amino acid > ] 
          [ -o < out pdb file > ]
================================================================================
r:  residue ID to be mutated, either residue ID or sequence position is must for 
    the operation to complete.
s:  sequence position to be mutated
a:  residue to mutate to (default mutation is ALA (A) )
o:  output file name ( default: <input>_mutate.pdb )
" 1>&2; }

cleanup(){ rm -f $tfile; }

halt(){ usage; cleanup; exit 1; }

resMut="A"

while getopts ":p:r:s:a:o:" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
        ;;
     o)
        outfile="$(readlink -f ${OPTARG})";
        ;;
     r)
        resId="${OPTARG}";
        ;;
     s)
        seqId="${OPTARG}";
        ;;
     a)
        resMut="${OPTARG}";
        ;;
     *)
        usage;
        exit 0;
        ;;
 esac
done
shift $((OPTIND-1));


if [ -z ${pdbfile} ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ] || [ $( $findNbr -f "$pdbfile" -C | wc -l) -ne 1 ] || [ "$(echo $pdbfile | grep -E ".pdb$")" != "$pdbfile" ]
then
   halt;
fi

if [ -z $outfile ]
then 
  outfile="${pdbfile%.pdb}_mutant.pdb"
fi

if [ -z ${outfile} ] || [ -d "${outfile}" ]
then
   halt;
fi

if [ "$(echo "ACDEFGHIKLMNPQRSTVWY" | grep $resMut)" == "" ]
then
   halt;
fi

if [ -z $resId ] && [ -z $seqId ] 
then
   halt;
fi

tfile=$(tempfile)

if [ ! -z $resId ]
then
  var="resid=$resId"
else
  var="seqid=$seqId"
fi

awk -v$var -vmut=$resMut 'BEGIN{ 
         resMap["ALA"] = "A";
         resMap["CYS"] = "C";
         resMap["ASP"] = "D";
         resMap["GLU"] = "E";
         resMap["PHE"] = "F";
         resMap["GLY"] = "G";
         resMap["HIS"] = "H";
         resMap["ILE"] = "I";
         resMap["LYS"] = "K";
         resMap["LEU"] = "L";
         resMap["MET"] = "M";
         resMap["ASN"] = "N";
         resMap["PRO"] = "P";
         resMap["GLN"] = "Q";
         resMap["ARG"] = "R";
         resMap["SER"] = "S";
         resMap["THR"] = "T";
         resMap["VAL"] = "V";
         resMap["TRP"] = "W";
         resMap["TYR"] = "Y";
         counter=1; 
     } 
     /^ATOM/ { 
        resId = sprintf( "%d", substr($0,23,4) ); 
        resName = substr($0,18,3); 
        if(residue[resId] == 0){ 
          residue[resId] = resMap[resName]; 
          rescount[counter] = resId;
          counter++;
        } 
      }
      END{
        for(i=1; i < counter; ++i){
          aa = residue[ rescount[i] ];
          if( resid == rescount[i] || seqid == i ){
             aa = mut; 
          }
          printf("%s", aa );
        }
        printf( "\n" );
      }' "$pdbfile" > $tfile


${scwrl} -s "${tfile}" -i "${pdbfile}" -o "${outfile}" 2>/dev/null 1>&2


cleanup;



