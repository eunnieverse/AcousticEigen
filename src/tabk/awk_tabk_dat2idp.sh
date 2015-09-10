#!/bin/bash 
########################################################################
### AcousticEigen project
### Yoonkyung Eunnie Lee, 2015.09.10
### convert tabk file from DAT to IDP form 
########################################################################
### initial format: 
### 67    1.128510    579829.671670    692245.615160
### changed format: 
### tabk(67)=579829.671670+(692245.615160i);    // Vec_kaLpi=1.128510
########################################################################
echo -n "Enter the name of the dat file and press [ENTER]: "
read DATFILENAME
FILEBASE=${DATFILENAME%.dat}

awk 'BEGIN {i=0;} { printf "tabk(%d)=%f+(%fi);    // Vec_kaLpi=%f\n",$1,$3,$4,$2;} {i=i+1;}' ${FILEBASE}.dat >> ${FILEBASE}.idp
