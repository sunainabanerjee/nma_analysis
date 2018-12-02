#! /bin/bash

if [ $# -lt 1 ] || [ ! -f $1 ]
then
  echo "Usage: $0 <mode files>" 1>&2;
  exit 1;
fi

pdbs=""

for pdb in $*
do
  if [ -f $pdb ]
  then
    pdbs="$pdbs $pdb"
  fi
done

awk '/^ATOM/ && substr($0,14,3) ~ /CA/ { 
          resId=sprintf("%d",substr($0,24,4)); 
          x=sprintf("%.3f",substr($0,30,8)); 
          y=sprintf("%.3f",substr($0,39,8)); 
          z=sprintf("%.3f",substr($0,47,8)); 
          if(X[resId] == ""){ 
             X[resId]=x; 
             Y[resId]=y; 
             Z[resId]=z; 
          } 
          RMS[resId] += (X[resId]-x)**2 + (Y[resId]-y)**2 + (Z[resId]-z)**2; 
          N[resId]++;
    }END{ 
      for(r in RMS) 
          printf("%d,%f\n",r,sqrt(RMS[r]/N[r]));
    }' $pdbs | sort -t, -k1n


