#!/bin/bash 
########################################################################
### AcousticEigen project
### Yoonkyung Eunnie Lee, 2015.09.10
### convert tabk file from IDP to DAT form 
########################################################################
### initial format: 
### tabk(67)=579829.671670+(692245.615160i);    // Vec_kaLpi=1.128510
### changed format: 
### 67    1.128510    579829.671670    692245.615160
########################################################################
echo -n "Enter the name of the idp file and press [ENTER]: "
read IDPFILENAME

FILEBASE=${IDPFILENAME%.idp}

cp  ${FILEBASE}.idp  ${FILEBASE}.txt 

sed -i 's/tabk(//g'  ${FILEBASE}.txt
sed -i 's/)=/    /g' ${FILEBASE}.txt
sed -i 's/+(/    /g' ${FILEBASE}.txt
sed -i 's/i)\;    \/\/ Vec_kaLpi=/    /g' ${FILEBASE}.txt

awk 'BEGIN {i=0;} { printf "%d    %f    %f    %f\n",$1,$4,$2,$3;} {i=i+1;}' ${FILEBASE}.txt >> ${FILEBASE}.dat

rm ${FILEBASE}.txt

