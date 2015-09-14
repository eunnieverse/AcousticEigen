%---------------------------------------------------------------------
%- plot tabk, main function 
%- Yoonkyung Eunnie Lee 2015. 09. 14
%---------------------------------------------------------------------
clear all; close all; 
addpath('/home2/eunnie12/Work/AcousticEigen_Run/matlab/');
addpath('/home2/eunnie12/Work/AcousticEigen_Run/matlab/utils');

filelist = {'NewtCyl_Mode1_phi90_tabk.dat',...
            'NewtCyl_Mode2_phi90_tabk.dat',...
            'NewtCyl_Mode1_phi90_err002_final_k.dat',...            
            'NewtCyl_Mode2_phi90_err002_final_k.dat'
            }; 

savejpg = 1; 
saveeps = 1; 

for ii=1:length(filelist)
    fignum = 1; 
    filename = filelist{ii}; 
    filebase = filename(1:length(filename)-4);
    cfig = plotk(filename);   
    plotsave(cfig,filebase,fignum,savejpg,saveeps); 
end