#! /bin/bash

exedir=$(dirname $(readlink -f $0))
consvExe="${exedir}/conservation_from_sequence_alignment.sh"

if [ ! -x "${consvExe}" ]
then
  echo "Error: installation error" 1>&2;
  exit 1;
fi 


cleanup(){
  rm -f $tfile $rfile $xfile
}

halt(){
  cleanup;
  echo "Error: in execution" 1>&2;
  exit 1;
}

usage(){
 echo "
------------------------------------------------------------------
Usage: $0 -a <sequence alignment file> 
          -t <sequence tag> 
          [-i <conservation cutoff>]
          [-m <minimum number of substitution>]
          [-c <normalize mututal information cutoff>]
------------------------------------------------------------------
" 1>&2;
}

min_substitution=3
consv_cutoff=60.
nmi_cutoff=0.8

while getopts ":a:t:i:m:c:" o;
do
  case "${o}" in
     a)
        alignment_file=$(readlink -f ${OPTARG})
        ;;
     t)
        tag=${OPTARG};
        ;;
     i)
        consv_cutoff=${OPTARG};
        ;;
     c)
        nmi_cutoff=${OPTARG};
        ;;
     m)
        min_substitution=${OPTARG};
        ;;
     *)
        usage;
        exit 1;
        ;;
  esac
done
shift $((OPTIND-1))


if [ -z $alignment_file ] || [ ! -f "${alignment_file}" ]
then
  usage;
  exit 1;
fi

if [ -z $tag ] || [ $(grep -E "^>.*${tag}" ${alignment_file} | wc -l) -ne 1 ]
then
  usage;
  exit 1;
fi

if [ -z $consv_cutoff ] || [ "$(echo $consv_cutoff | grep -P '^\d+\.?\d*$')" != "${consv_cutoff}" ]
then
  usage;
  exit 1;
fi

if [ -z $min_substitution ] || [ "$(echo $min_substitution | grep -P '^\d+$')" != "${min_substitution}" ]
then
  usage;
  exit 1;
fi

tfile=$(tempfile)
rfile=$(tempfile)
xfile=$(tempfile)

awk '/^>/ && NR > 1 {print ""} !/^>/{printf "%s", $1}' $alignment_file  > $tfile

if [ $(awk '{print length($1)}' $tfile | sort -nu | wc -l) -ne 1 ]; then halt;fi 

seqLen=$(awk '{print length($1)}' $tfile | sort -nu)

pos=1;
while [ $pos -le $seqLen ]
do
   if [ $(awk -vpos=$pos 'NF==1 { x = substr($1,pos,1); if( x != "-" ) print x; } ' $tfile | sort -u | wc -l) -ge $min_substitution ]
   then
     echo $pos
   fi 
   pos=$((pos+1))
done > $xfile

"${consvExe}" -t "${tag}" -a "${alignment_file}" -c | awk -F, -vfile=$xfile -vconsvCutoff=${consv_cutoff} 'BEGIN{
                                                               while(getline line < file)
                                                                 allowed[line] = 1;
                                                            }
                                                            NF==3 && allowed[$1] == 1 && $3 < consvCutoff {
                                                              printf("%d,%.3f\n",$2,$3);
                                                            }' > ${rfile}

if [ $? -ne 0 ];then  halt; fi

pos=1

awk 'function mycmp(ia, va, ib, vb, sa, sb){
       if( int(ia) < int(ib) )    return -1;
       else return 1;
     }
    NF==1{ 
        pos[$1]= 1;
    } END{
       PROCINFO["sorted_in"] = "mycmp";
       for(i in pos)
         for(j in pos) 
            if( i > j) 
                printf("%d,%d\n",i,j);
    }' $rfile | while read line; do
   pos1=$(echo $line | awk -F, '{print $1}');
   pos2=$(echo $line | awk -F, '{print $2}');
   score=$(awk -vpos1=$pos1 -vpos2=$pos2 '{ 
            x = substr($1,pos1,1); 
            y = substr($1,pos2,1);
            
            if( x != "-" && y != "-"){
              xy= sprintf("%s%s",x,y); 
              p[x]++; 
              p[y]++;
              p1[x]++;
              p2[y]++;
              pxy[xy]++;
            } 
          }END{ 
            t = 0; 
            for(x in p) 
              t+= p[x]; 
            for(x in p) 
              p[x] = p[x] * 1.0/ t; 

            t = 0;
            for(x in p1)
              t += p1[x];
            for(x in p1)
              p1[x] = p1[x] * 1.0/t;
           
            t=0;
            for(x in p2)
              t += p2[x];
            for(x in p2)
              p2[x] = p2[x] * 1.0/t;

            t=0; 
            for(x in p) 
              for(y in p){ 
                xy=sprintf("%s%s",x,y); 
                t += pxy[xy];
              } 
            for(x in p)
               for( y in p){
                 xy = sprintf("%s%s", x, y);
                 pxy[xy] = pxy[xy]* 1.0/t;
               }
            score = 0;

            for(x in p)
             for(y in p){
               xy = sprintf("%s%s",x,y);
               if(pxy[xy] > 0.0){
                 score += pxy[xy] * log( pxy[xy]/(p1[x]*p2[y]) );
               }
             }
           
            h1 = 0;
            h2 = 0; 
            for(x in p1)
               h1 += -1. * p1[x] * log(p1[x]);
            for(x in p2)
               h2 += -1. * p2[x] * log(p2[x]);
           
            h = (h1 > h2)?h2:h1;

            printf( "%.4f\n", score/h );

          }' $tfile 2>/dev/null);
   if [ "$score" != "" ]
   then
      echo "$pos1,$pos2,$score";
   fi
done | awk -F, -vnmi_cutoff=${nmi_cutoff} '$3 > nmi_cutoff'

cleanup;
exit 0;
