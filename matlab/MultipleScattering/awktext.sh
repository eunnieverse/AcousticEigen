awk 'BEGIN {i=0;} { printf "tabk(%d)=%f+(%fi);    // Vec_kaLpi=%f\n",i,$2,$3,$1;} {i=i+1;}' tabk_to_ms_k.txt >tabk_to_ms_k_edp.txt
