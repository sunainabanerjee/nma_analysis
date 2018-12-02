#! /bin/bash

usage() {
  echo "Usage: $0 <list of one letter amino codes>" 1>&2;
}

if [ $# -eq 0 ]  
then
  usage;
  exit 1;
fi

for i in $*
do
  if [ "$(echo "$i" | grep -iP "^[A-Z]$")" !=  "$i" ]
  then
    echo "Error: invalid input: $i" 1>&2;
    usage;
    exit 1;
  fi
done

while [[ $# -gt 0 ]]
do
 echo $1 | awk 'BEGIN{ 
                 map["A"] = "ALA";
                 map["C"] = "CYS";
                 map["D"] = "ASP";
                 map["E"] = "GLU";
                 map["F"] = "PHE";
                 map["G"] = "GLY";
                 map["H"] = "HIS";
                 map["I"] = "ILE";
                 map["K"] = "LYS";
                 map["L"] = "LEU";
                 map["M"] = "MET";
                 map["N"] = "ASN";
                 map["P"] = "PRO";
                 map["Q"] = "GLN";
                 map["R"] = "ARG";
                 map["S"] = "SER";
                 map["T"] = "THR";
                 map["V"] = "VAL";
                 map["W"] = "TRP";
                 map["Y"] = "TYR";
               }
               {
                 for(i=1; i <= NF; ++i){
                   inp = toupper($i);
                   if( toupper(inp) ~ /^[A-Z]$/ &&  map[inp] != "" ){
                     printf("%s\n", map[inp]);
                   }else{
                     printf("%s\n", "UNK");
                   }
                     
                 }
               }'
  shift;
done
