#!/bin/bash

for dirname in Re_vxjr Re_vyjr Im_vxjr Im_vyjr Re_vxjc Re_vyjc Im_vxjc Im_vyjc Mag_vxjr Mag_vyjr Mag_vxjc Mag_vyjc Mag_vjr Mag_vjc Mag_Pjr Mag_Phc Mag_Tjr Mag_Tjc; 
do
  mkdir ${dirname}; 
  mv ${dirname}*.eps ${dirname}/.
done;
