%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
%%%%%%% Last updated on 2015-06-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set output plot conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extension = '.png'; %%default plot extension 
showplot = 1; 
    savefig = 1; 
    saveeps = 1; 
savetxt = 0; 
savemat = 0; 
prec = '%e'; %% precision 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set input directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basedir = pwd(); %% directory containing .m file 
basedir = '/home/eunnie12/Work/AcousticEigen_Run';
cd(basedir); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=293.; 			   % [K]
c0=331.+0.6*(T0-273.); % [m/s]
L=1e-5;                % [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load reference data (multiple scattering) 
load('MultipleScattering/ms.mat'); %c1,c2,k1,k2,klpi_ms
load('freefem/Local_correct/local.mat'); %cL, kL, klpi_l; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd('freefem/Mode1/');
%dirlist=dir('NewtCyl_Mode1_phi90_err005.txt');     %%for the first .txt file in the folder
cd('freefem/Mode1/');
dirlist=dir('NewtCyl_Mode1_phi90_err005.txt');     %%for the first .txt file in the folder

filename=dirlist.name;
fnsplit = strsplit(filename,'.'); %% fn= 1x2 cell, 'abcd_efgh' 'edp' 
fnext = char(fnsplit(length(fnsplit))); %% 'edp' 
filebase = char(fnsplit(1));

fid = fopen(filename);
   formatSpec = '%f %f %f %f %f';
   data = textscan(fid,formatSpec,'HeaderLines',0,'CollectOutput',1,'EmptyValue',0);
   data = data{1,1};
   fclose(fid); 
   %%% store different versions of x axis 
   klpi=data(:,1); %%% Omega is in the order of 2*pi*10^(-6)*/lambda.
   xlength=length(klpi);
   klpi_100=logspace(log10(0.05),log10(5),100);           
   Omega = klpi*c0*pi/L;
   Omega_100=klpi_100*c0*pi/L; 
             
   Rek = data(:,   2);
   Imk = data(:,   3);
   Rec = data(:,   4);
   Imc = data(:,   5);
   clear data; 
   
%k1 = Omega./c1; %%% K1 in .mat file was wrong. 
%kL = Omega./cL; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(showplot==1)
  xlab='kL/\pi, dimensionless frequency';
  xlimit1=[min(klpi) 5]; %%limit is 0.05~5 
  xtick1=[0.05; 0.25; 0.5; 2.5; 5];
  xlimit2=[min(klpi) 2]; %%limit is 0.05~5 
  xtick2=[0.05; 0.25; 0.5; 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig=figure(); %%%Re[k]
    semilogx(klpi,Rek,'ko','LineWidth',2);
    hold on; semilogx(klpi_ms,real(k1),'b',klpi_ms,real(k2),'g--','LineWidth',2); 
    hold on; semilogx(klpi_ms,real(kL),'r','LineWidth',2); 
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    ylabel('wavevector k [1/m]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
    savefigname=strcat('Rek_5_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

cfig=figure(); %%%Im[k]
    semilogx(klpi,Imk,'ko','LineWidth',2);
    hold on; semilogx(klpi_ms,imag(k1),'b',klpi_ms,imag(k2),'g--','LineWidth',2); 
    hold on; semilogx(klpi_ms,imag(kL),'r','LineWidth',2); 
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    ylabel('wavevector k [1/m]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
    savefigname=strcat('Imk_5_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
    
cfig=figure(); %%%Re[k], xlimit2
    semilogx(klpi,Rek,'ko','LineWidth',2);
    hold on; semilogx(klpi_ms,real(k1),'b',klpi_ms,real(k2),'g--','LineWidth',2); 
    hold on; semilogx(klpi_ms,real(kL),'r','LineWidth',2); 
    set(gca,'XLim',xlimit2,'XTick',xtick2); xlabel(xlab);
    ylabel('wavevector k [1/m]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
    savefigname=strcat('Rek_2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

cfig=figure(); %%%Im[k], xlimit2
    semilogx(klpi,Imk,'ko','LineWidth',2);
    hold on; semilogx(klpi_ms,imag(k1),'b',klpi_ms,imag(k2),'g--','LineWidth',2); 
    hold on; semilogx(klpi_ms,imag(kL),'r','LineWidth',2); 
    set(gca,'XLim',xlimit2,'XTick',xtick2); xlabel(xlab);
    ylabel('wavevector k [1/m]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
    savefigname=strcat('Imk_2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x range from 0.05 ~ 5.0 
cfig = figure();
        semilogx(klpi,Rec,'ko','LineWidth',2);
        hold on; semilogx(klpi_ms,real(c1),'b',klpi_ms,real(c2),'g--','LineWidth',2); 
        hold on; semilogx(klpi_ms,real(cL),'r','LineWidth',2); 
        set(gca,'XLim',xlimit1,'XTick',xtick1);
        xlabel(xlab);
        ylabel('phase velocity [m/s]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
        savefigname=strcat('Rec_5_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

    cfig = figure();
        semilogx(klpi,Imc,'ko','LineWidth',2);
        hold on; semilogx(klpi_ms,imag(c1),'b',klpi_ms,imag(c2),'g--','LineWidth',2); 
        hold on; semilogx(klpi_ms,imag(cL),'r','LineWidth',2); 
        set(gca,'XLim',xlimit1,'XTick',xtick1);
        xlabel(xlab);
        ylabel('phase velocity [m/s]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
        savefigname=strcat('Imc_5_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x range from 0.05 ~ 2.0 
cfig = figure();
        semilogx(klpi,Rec,'ko','LineWidth',2);
        hold on; semilogx(klpi_ms,real(c1),'b',klpi_ms,real(c2),'g--','LineWidth',2); 
        hold on; semilogx(klpi_ms,real(cL),'r','LineWidth',2); 
        set(gca,'XLim',xlimit2,'XTick',xtick2);
        xlabel(xlab);
        ylabel('phase velocity [m/s]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
        savefigname=strcat('Rec_2_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

    cfig = figure();
        semilogx(klpi,Imc,'ko','LineWidth',2);
        hold on; semilogx(klpi_ms,imag(c1),'b',klpi_ms,imag(c2),'g--','LineWidth',2); 
        hold on; semilogx(klpi_ms,imag(cL),'r','LineWidth',2); 
        set(gca,'XLim',xlimit2,'XTick',xtick2);
        xlabel(xlab);
        ylabel('phase velocity [m/s]');
    h=legend('FEM','nonlocal1','nonlocal2','local'); set(h,'Location','northeastoutside');
        savefigname=strcat('Imc_2_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
end %%if(showplot)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Save as separate TXT File if asked 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(savetxt==1)
%    outname= strcat('tabk','_',filebase,'.txt'); 
%    dlmwrite(outname,[klpi',Rek,Imk],'delimiter','\t','precision',prec);
% end %%if(savetxt) 
  cd (basedir); 