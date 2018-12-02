#! /bin/bash

usage() {
  echo "Usage: $0 <list of 3 letter amino codes>" 1>&2;
}

if [ $# -eq 0 ]  
then
  usage;
  exit 1;
fi

for i in $*
do
  if [ "$(echo "$i" | grep -iP "^[A-Z]{3}$")" !=  "$i" ]
  then
    echo "Error: invalid input: $i" 1>&2;
    usage;
    exit 1;
  fi
done

while [[ $# -gt 0 ]]
do
 echo $1 | awk 'BEGIN{ 
                 map["ALA"] = "A";
                 map["CYS"] = "C";
                 map["ASP"] = "D";
                 map["GLU"] = "E";
                 map["PHE"] = "F";
                 map["GLY"] = "G";
                 map["HIS"] = "H";
                 map["ILE"] = "I";
                 map["LYS"] = "K";
                 map["LEU"] = "L";
                 map["MET"] = "M";
                 map["ASN"] = "N";
                 map["PRO"] = "P";
                 map["GLN"] = "Q";
                 map["ARG"] = "R";
                 map["SER"] = "S";
                 map["THR"] = "T";
                 map["VAL"] = "V";
                 map["TRP"] = "W";
                 map["TYR"] = "Y";
               }
               {
                 for(i=1; i <= NF; ++i){
                   inp = toupper($i);
                   if( toupper(inp) ~ /^[A-Z]{3}$/ &&  map[inp] != "" ){
                     printf("%s\n", map[inp]);
                   }else{
                     printf("%s\n", "X");
                   }
                     
                 }
               }'
  shift;
done
