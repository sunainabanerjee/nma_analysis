#! /bin/bash


usage(){ echo "Usage: $0 -p <pdb file> [-c <chain ID>] [-r <distance cutoff>] [-n <minimum atom contacts>]" 1>&2;}


while getopts ":p:c:r:n:" o;
do
  case "${o}" in
     p)
        pdbfile=$(readlink -f ${OPTARG})
        ;;
     c)
        chain=${OPTARG};
        ;;
     r)
        radius=${OPTARG};
        ;;
     n)
        cutoff=${OPTARG};
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${pdbfile}" ] || [ ! -f "${pdbfile}" ] || [ $(awk '/^ATOM/' "${pdbfile}"| wc -l) -lt 4 ]
then
  usage;
  exit 1;
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
  usage;
  exit 1;
fi

if [ -z "${radius}" ] || [ "$(echo $radius | grep -P '\d+\.?\d*')" != "${radius}" ]
then
  radius=4.5
fi

if [ -z "${cutoff}" ] || [ "$(echo $cutoff | grep -P '\d+')" != "${radius}" ]
then
  cutoff=3
fi

xfile=$(tempfile)

if [ $(find_neighbor_pdb -f "${pdbfile}" -C 2>/dev/null | awk -F: -vchain=$chain '$1 == chain' | wc -l)  != 1 ]
then
  echo "Error: chain ($chain) not available in the pdb [$pdbfile]" 1>&2;
  exit 1;
fi

find_neighbor_pdb -f "${pdbfile}" -c "${chain}" -e 2>/dev/null 1>"${xfile}"

awk 'substr($0,14,3) ~ "CA" && /^ATOM / { 
        resId=substr($0,23,4); 
        gsub(" ", "",resId); 
        resName=substr($0,18,3); 
        printf("%s%d\n", resName,resId); 
     }' ${xfile} | 
            uniq | while read resId; 
                   do 
                   find_neighbor_pdb -f "${xfile}" -r $radius -S ${resId}:${chain} 2>/dev/null | 
                             awk -vtarget=$resId -vminContact=$cutoff '/^ATOM / { 
                                                                        targetRes=substr(target,4);
									resId=substr($0,23,4); 
                                                    	                gsub(" ", "",resId); 
                                                                        resName=substr($0,18,3); 
                                                                        diff = ( resId > targetRes )?(resId - targetRes):(targetRes - resId);
                                                                        if( diff > 1 ){
                                                                         contact=resName""resId; 
                                                                         ccount[contact]++;
                                                                        }
                                                                  }
                                                                  END{
                                                                    for( c in ccount)
                                                                      if( ccount[c] >= minContact )
                                                                        printf("%s,%s\n",target,c);
                                                                  }' ; 
                   done

rm -f ${xfile}
