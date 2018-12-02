#! /bin/bash

usage(){
  echo "Usage: $0 -1 <pdb file1> -2 <pdb file2> -a <pdb file1 chain id> -b <pdb file-2 chain id>" 1>&2;
}

while getopts ":1:2:a:b:" o;
do
  case "${o}" in
    1) pdbfile1=$(readlink -f ${OPTARG});
       ;;
    2) pdbfile2=$(readlink -f ${OPTARG});
       ;;
    a) chain1=${OPTARG};
       ;;
    b) chain2=${OPTARG};
       ;;
    *) usage;
       exit 1;
  esac
done
shift $((OPTIND-1))

if [ -z "${pdbfile1}" ] || [ ! -f "${pdbfile1}" ] || [ $(awk '/^ATOM/' "${pdbfile1}"| wc -l) -lt 4 ]
then
  usage;
  exit 1
fi 

if [ -z "${pdbfile2}" ] || [ ! -f "${pdbfile2}" ] || [ $(awk '/^ATOM/' "${pdbfile2}"| wc -l) -lt 4 ]
then
  usage;
  exit 1
fi 

if [ -z "${chain1}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile1}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain1=$(find_neighbor_pdb -f "${pdbfile1}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi

if [ -z "${chain2}" ]
then
   if [ $( find_neighbor_pdb -f "${pdbfile2}" -C 2>/dev/null | wc -l) -eq 1 ]
   then
     chain2=$(find_neighbor_pdb -f "${pdbfile2}" -C 2>/dev/null | awk -F: '{print $1}')
   fi
fi

if [ $(echo -n "$chain1" | wc -c) -ne 1 ] || [ $(echo -n "$chain2" | wc -c) -ne 1 ]
then
  usage;
  exit 1;
fi

tdir="/tmp/salign_$$"

if [ -d "$tdir" ]
then
  rm -fr "${tdir}" || exit 1
fi

mkdir -p "${tdir}"

pfile1="${tdir}/pdbfile1.pdb"
pfile2="${tdir}/pdbfile2.pdb"

find_neighbor_pdb -f "${pdbfile1}" -c "${chain1}" -e 2>/dev/null > "${pfile1}"
find_neighbor_pdb -f "${pdbfile2}" -c "${chain2}" -e 2>/dev/null > "${pfile2}"

currDir=$(pwd)
cd $tdir

Matt -o map "${pfile1}" "${pfile2}"

alnfasta="map.fasta"

if [ -f "${alnfasta}" ]
then
  awk -vhdr1="pdbfile1" -vhdr2="pdbfile2" '/^>/{ 
                                              key=substr($1,2); 
                                              split(key,arr,":"); 
                                              hdr=arr[1]; 
                                           } !/^>/{ 
                                              map[hdr]=map[hdr]""$1;
                                           } END{ 
                                              n=length(map[hdr1]); 
                                              for(i=1; i<=n; ++i){ 
                                               c1=substr(map[hdr1],i,1); 
                                               c2=substr(map[hdr2],i,1); 
                                               if(c1 != "-") p1++; 
                                               if(c2 != "-") p2++; 
                                               if(c1 != "-" && c2 != "-") 
                                                  printf("%d,%d,%d\n",i,p1,p2); 
                                              }  
                                          }' map.fasta
fi

cd ${currDir}
rm -fr "${tdir}"
