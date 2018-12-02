#! /bin/bash

exedir=$(dirname $(readlink -f $0));

modelBurstExe="${exedir}/burst_model.sh"

for e in "${modelBurstExe}" 
do
  if [ ! -x "${e}" ]
  then
    echo "Error: setup error [missing ${e}]" 1>&2;
    exit 1;
  fi
done


usage(){
  echo "
===================================================
Usage: $0 
          -p <pdb model file>
          -l <residue pair list>
         [-c <pdb chain>]
         [-s]
===================================================
p: pdb model must contain multiple model of same protein
l: list of residue pairs. The residue pairs must be specified
   in the following format, 
c: in case of models with multiple chain, provide single letter
   chain identifier.
s: switch on scrwl
>>>EXAMPLE
ARG262,LEU250
LEU289,ALA112
>>>END
" 1>&2;
}


cleanup(){ rm -fr $tfile $tdir; }

halt(){ cleanup; usage; exit 1; }

rotamer=""

extension='mode'
while getopts ":p:l:c:s" o;
do
  case "${o}" in
     p)
        pdbfile="$(readlink -f ${OPTARG})";
        ;;
     c)
        chain=${OPTARG};
        ;;
     l)
        infile="$(readlink -f ${OPTARG})";
        ;;
     s)
        rotamer="-s";
        ;;
     *)
        halt;
        ;;
  esac
done
shift $((OPTIND-1))

tfile=$(tempfile)
tdir="/tmp/pair_fluctuation_$$"


if [ -z "${pdbfile}" ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ] || [ $(grep -P '^MODEL\s+\d+' "${pdbfile}" | wc -l ) -lt 3 ]
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

if [ -z $chain ] || [ $(echo -n "$chain" | wc -c) -ne 1 ]
then
  halt;
fi

if [ -z $infile ] || [ ! -f "${infile}" ] || [ "$(awk -F, '$1 ~ /^[A-Z]{3}-?[0-9]+$/ && $2 ~  /^[A-Z]{3}-?[0-9]+$/' "${infile}" | wc -l)"  !=  "$( cat "${infile}" | wc -l)" ]
then
  halt;
fi


if [ ! -d "$tdir" ]
then
  mkdir -p "$tdir" || halt;
fi


"${modelBurstExe}" -p "${pdbfile}" -o "${tdir}" -e "${extension}" "${rotamer}"

currdir=$(pwd)
cd "${tdir}"

pdbs=$(ls ${extension}*pdb | sort -t_ -k2n)

awk -vfile="${infile}" -vchain=$chain -vext=${extension} 'BEGIN{
       while( getline line < file ){
         split( line, arr, ",");
         if( length(arr) >= 2 ){
            node[arr[1]] = 1;
            node[arr[2]] = 1;
         }
       }
     }
     /^ATOM / && substr($0,22,1) == chain {
       tag=FILENAME;
       gsub(ext, "", tag);
       gsub("_", "", tag);
       gsub(".pdb", "", tag);
       resId=sprintf("%d", substr($0,23,4));
       resName=substr($0,18,3); 
       restag=sprintf( "%s%d", resName, resId); 
       atomName=substr($0,13,4);
       x = sprintf("%.3f", substr($0,31,8));
       y = sprintf("%.3f", substr($0,39,8));
       z = sprintf("%.3f", substr($0,47,8));
       if( node[restag] == 1 )
         printf("%2d %3s %4d %4s %8.3f %8.3f %8.3f\n", tag,resName,resId,atomName,x,y,z);
     }'  $pdbs > $tfile
cd $currdir

awk -vfile="${infile}" 'BEGIN{ 
      while(getline line < file){ 
         split(line,arr,","); 
         if( length(arr) >=2 ){
           edge[arr[1],arr[2]]=1;
         }   
      } 
     } NF==7{ 
       restag=sprintf("%s%d", $2,$3); 
       mode[$1]=1; 
       residue[restag]=1; 
       atom[$4]=1; 
       resatom[restag,$4]=1; 
       X[$1,restag,$4]=$5; 
       Y[$1,restag,$4]=$6; 
       Z[$1,restag,$4]=$7;
     }END{ 
        for( r1 in residue ){ 
          for(r2 in residue ){ 
            if(edge[r1,r2] == 1){ 
              x = 0; x2= 0; nx = 0;
              for(m in mode){
                d   = sqrt((X[m,r1,"CA"] - X[m,r2,"CA"])**2 + (Y[m,r1,"CA"] - Y[m,r2,"CA"])**2 + (Z[m,r1,"CA"] - Z[m,r2,"CA"])**2);
                x  += d;   x2 += d*d;   nx++;
              }

              for( m in mode){
                for(a1 in atom){
                  if( resatom[r1,a1] == 1){
                    for(a2 in atom){
                      if( resatom[r2,a2] == 1 ){
                         d = sqrt((X[m,r1,a1] - X[m,r2,a2])**2 + (Y[m,r1,a1] - Y[m,r2,a2])**2 + (Z[m,r1,a1] - Z[m,r2,a2])**2);
                         dz[a1,a2] += d;  dz2[a1,a2] += d*d;  nz[a1,a2]++;
                      }
                    }
                  }
                }
              }

              y=0;
              for(a1 in atom){
                if( resatom[r1,a1] == 1 ){
                  for( a2 in atom){
                     if(resatom[r2,a2] == 1 ){
                        yk = dz[a1,a2]/nz[a1,a2];
                        if( yk > y ){
                          na1 = a1;  na2 = a2; y = yk;    
                        }
                     }
                  }
                }
              }

              z = 100;
              for(a1 in atom){
                if( resatom[r1,a1] == 1 && a1 ~ /N|O|S/ ){
                  for( a2 in atom){
                     if(resatom[r2,a2] == 1 && a2 ~ /N|O|S/ ){
                        zk = dz[a1,a2]/nz[a1,a2];
                        if( zk < z ){
                          ma1 = a1;  ma2 = a2; z = zk;    
                        }
                     }
                  }
                }
              }
              
              printf("%s,%s,%.3f,%.3f,%s,%s,%.3f,%.3f,%s,%s,%.3f,%.3f\n",r1,r2,
                                                       x/nx,
                                                       sqrt((x2/nx) - (x/nx)**2),
                                                       na1,na2,y,
                                                       sqrt((dz2[na1,na2]/nz[na1,na2]) - y**2),
                                                       ma1,ma2,z,
                                                       sqrt( (dz2[ma1,ma2]/nz[ma1,ma2]) - z**2));
              delete dz;
              delete dz2;
              delete nz;
            }   
          }  
        } 
     }' $tfile

cleanup;
exit 0;
