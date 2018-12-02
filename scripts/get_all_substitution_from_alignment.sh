#! /bin/bash

exedir=$(dirname $(readlink -f $0))
consvExe="${exedir}/conservation_from_sequence_alignment.sh"

if [ ! -x "${consvExe}" ]
then
  echo "Error: installation error" 1>&2;
  exit 1;
fi

cleanup(){ rm -f $xfile $tfile; }

usage(){ echo "Usage: $0 -t <tag name> -a <alignment file>" 1>&2; }

halt() { cleanup; exit 1;}

while getopts ":t:a:c" o;
do
  case "${o}" in
     t)
        tag="${OPTARG}"
        ;;
     a)
        alignment_file=$(readlink -f "${OPTARG}");
        ;;
     *)
        usage;
        exit 1;
        ;;
  esac
done
shift $((OPTIND-1))

tfile=$(tempfile)
xfile=$(tempfile)

if [ -z ${alignment_file} ] || [ ! -f "${alignment_file=}" ]; then usage; halt; fi

if [ -z ${tag} ] || [ $(grep -P "^>.*${tag}" "${alignment_file}"| wc -l) -ne 1 ]; then usgae; halt; fi

"${consvExe}" -a "${alignment_file}" -t $tag | awk -F, '{print $1}' > $xfile

awk '/^>/ && NR > 1 {print ""} !/^>/{printf "%s", $1}' "${alignment_file}"  > $tfile

if [ $(awk '{print length($1)}' $tfile | sort -nu | wc -l) -ne 1 ]
then
  echo "Error: improper alignment file ${alignment_file}" 1>&2;
  exit 1;
fi

seqLen=$(awk '{print length($1)}' $tfile | sort -nu)

pos=1

for pos in $(cat $xfile)
do
  echo ${pos} $(awk -vpos=$pos 'NF==1 && length($1) >= pos { x=substr($1,pos,1); if(x != "-") print x }' $tfile | sort -u);
  pos=$((pos+1))
done

cleanup
exit 0
