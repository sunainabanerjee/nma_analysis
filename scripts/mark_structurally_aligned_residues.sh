#! /bin/bash

exedir=$(dirname $(readlink -f $0))
mapexe="${exedir}/map_seqpos2restag.sh"


if [ ! -x "${mapexe}" ]
then
  echo "Error: installation error!" 1>&2;
  exit 1;
fi 


usage(){ echo "Usage: $0 <list of pdbs..>" 1>&2; }


if [ $# -lt 2 ]
then
  usage;
  exit 1; 
fi

counter=$#
basefile="/tmp/tfile_$$";
#tmpdir="/tmp/tdir_$$";
tmpdir="/tmp/tdir";
baseCounter=0
for i in $*; 
do
   if [ ! -f "$i" ] || [ $(awk '/^ATOM/' "$i" |wc -l) -lt 4 ] || [ $(find_neighbor_pdb -f "${i}" -C 2>/dev/null| wc -l ) -gt 1 ]
   then
      echo "Error: improper pdb ${i}" 1>&2;
      exit 1;
   fi
   tfile[$baseCounter]=$(tempfile)
   ifile[$baseCounter]="$(readlink -f "$i")"
   iname[$baseCounter]="$(basename "$i")"
   pname[$baseCounter]="$(basename "$i" | sed -e 's/.pdb$//')"
   $mapexe -p "${i}" > "${tfile[$baseCounter]}"
   baseCounter=$((baseCounter+1));
done



if [ ! -d "${tmpdir}" ]
then
  mkdir -p "${tmpdir}"
fi

currDir=$(pwd)
cd "${tmpdir}"

for i in ${ifile[*]};
do
  cp "${i}" .
done

Matt -o map ${iname[*]} >/dev/null 2>&1

alnfasta="map.fasta"

if [ -f "${alnfasta}" ]
then
  awk '/^>/{ 
            key=substr($1,2); 
            split(key,arr,":"); 
            hdr=arr[1]; 
       } !/^>/{ 
           map[hdr]=map[hdr]""$1;
       } END{ 
          n=length(map[hdr]); 
          counter = 0;
          for(j in map){
           pos[j] = 0;
           if( counter > 0) printf(",");
           printf("%s", j);
           counter++;
          }
          printf("\n");
          for(i=1; i<=n; ++i){
            p = 0;
            for( j in map ){
              h[j] = 0;
              c = substr(map[j], i, 1);
              if( c != "-" ){
                 p++;
                 h[j]++;
                 pos[j]++;
              }
            } 
            if( p > 1 ){
               counter = 0;
               for(j in map){
                 if(counter > 0) printf(",");
                 u = (h[j] > 0)?pos[j]:"-";
                 printf("%s",u);
                 counter++;
               }
               printf("\n");
            }
          }  
       }' map.fasta
fi

cd "${currDir}"
rm -fr "${tmpdir}"
rm -f ${tfile[*]}
