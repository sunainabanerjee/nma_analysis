#! /bin/bash


usage(){ echo "Usage: $0 -t <tag name> -a <alignment file> [-c]" 1>&2; }


consv=0
while getopts ":t:a:c" o;
do
  case "${o}" in
     t)
        tag="${OPTARG}"
        ;;
     a)
        alignment_file=$(readlink -f "${OPTARG}");
        ;;
     c)
        consv=1;
        ;;
     *)
        usage;
        exit 0;
        ;;
  esac
done
shift $((OPTIND-1))

if [ -z ${alignment_file} ] || [ ! -f "${alignment_file}" ]
then
   usage;
   exit 1;
fi

if [ -z ${tag} ] || [ $(grep -P "^>.*${tag}" "${alignment_file}"| wc -l) -ne 1 ]
then
  usage;
  exit 1;
fi

tfile=$(tempfile)
rfile=$(tempfile)
qfile=$(tempfile)
xfile=$(tempfile)
ufile=$(tempfile)
pfile=$(tempfile)


awk '/^>/ && NR > 1 {print ""} !/^>/{printf "%s", $1}' "${alignment_file}"  > $tfile

if [ $(awk '{print length($1)}' $tfile | sort -nu | wc -l) -ne 1 ]
then
  echo "Error: improper alignment file ${alignment_file}" 1>&2;
  exit 1;
fi

seqLen=$(awk '{print length($1)}' $tfile | sort -nu)

pos=1

while [ $pos -le $seqLen ]
do
	score=$(awk -vpos=${pos}  '{print substr($1,pos,1)}'  $tfile | 
                                                             uniq -c |
                                                             awk 'NF == 2 && $2 != "-" {
                                                                  pos[$2] = $1; 
                                                                  count++;
                                                              }END{ 
                                                                score = 0; 
                                                                if(count > 1){ 
                                                                  t=0; 
                                                                  for(i in pos) 
                      	                                             t+= pos[i];  
                                                                  n = length(pos);
                                                                  if( pos["-"] > 0 )
                                                                     n = n - 1;
                                                                  contrib = pos["-"] * (1.0/n);
	                                                          for(i in pos){
                                                                    if( i != "-"){ 
								      p = (pos[i] + contrib)/t; 
                                                                      score += -1. * p * log(p);
                                                                    }
                                                                  }
								 }
                                                                 if( n > 1 )
                                                                     printf("%.4f", score/log(n));
                                                                 else
                                                                     printf("%.4f", 0);
                          
								}');
       if [ "${score}" != "" ]
       then
            echo "$pos,$score";
       fi
       pos=$((pos+1))
done  > ${rfile}

awk -F, 'NF==2{print $1;}' "${rfile}"  |
                             sort -nu | while read pos; 
                                         do
                                            spos=$(sed -n "/>.*${tag}/,/>/p" "${alignment_file}" |
                                                                                    grep -v "^>" |
                                                                         awk '{printf "%s", $1}' |
                                                       awk -vpos=$pos 'substr($1,pos,1) != "-" {
                                                              printf "%s",substr($1,1,pos)
                                                       }' | sed -e 's/-//g' | wc -c);
                                            if [ ${spos} -ne 0 ];
                                            then
                                               echo "$pos,$spos";
                                            fi
                                            done > "${pfile}"

awk -vfile="${pfile}" -vconsv=${consv} -F, 'BEGIN{
                            while(getline line < file){ 
                               split(line, arr, ",");
                               if( length(arr) == 2 )
                                   pos[arr[1]] = arr[2]; 
                            }
                          }
                          NF==2 && pos[$1] != 0{
                             score = (consv == 1)?(1 - $2):$2;
                             score = score * 100;
                             printf("%d,%d,%.3f\n",$1,pos[$1],score);
                          }' "${rfile}"

rm -f $tfile $rfile $qfile $xfile $ufile $pfile
