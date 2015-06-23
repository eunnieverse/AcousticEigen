%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
%%%%%%% Last updated on 2015-06-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set output plot conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extension = '.png'; %%default plot extension 
showplot = 1; 
    savefig = 1; 
    saveeps = 0; 
savetxt = 1; 
savemat = 0; 
unitconst = 2*pi / 1240;
prec = '%e'; %% precision 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set input directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basedir = pwd(); %% directory containing .m file 
basedir = '/home/eunnie12/Work/AcousticEigen';
cd(basedir); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=293.; 			   % [K]
c0=331.+0.6*(T0-273.); % [m/s]
L=1e-5;                % [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load multiple scattering data from 1000-pt mat file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load reference data (multiple scattering) 
%%% 1st mode from *.mat
%%% 2nd mode from tabk_to_ms.txt 
load('matlab/Navid_Phi90_DATA_PTS100.mat'); 
%klpi_ms = Vec_kalpi; 
data =dlmread('matlab/tabk_to_ms_c.txt'); 
klpi_ms=data(:,1);
c2= data(:,2) + i*data(:,3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('freefem/Mode1/');
dirlist=dir('NewtCyl_Mode1_phi90_err002*.txt');     %%for the first .txt file in the folder
filename=dirlist.name
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
   klpi_316=logspace(log10(0.05),log10(5),316); 
   klpi_100=logspace(log10(0.05),log10(5),100);           
   Omega = klpi*c0*pi/L;
   Omega_316=klpi_316*c0*pi/L; 
   Omega_100=klpi_100*c0*pi/L; 
             
   Rek = data(:,   2);
   Imk = data(:,   3);
   Rec = data(:,   4);
   Imc = data(:,   5);
   clear data; 
            
% filename2= 'Valeur_k__Phi90_ka_complexe.txt'; 
%             fid = fopen(filename2);
%             formatSpec = '%f %f';
%             data = textscan(fid,formatSpec,'HeaderLines',0,'CollectOutput',1,'EmptyValue',0);
%             data = data{1,1};
%             fclose(fid); 
%             
%             Rek_ms=data(:,1);
%             Imk_ms=data(:,2);
%             k_ms=Rek_ms+Imk_ms*i;
%             c_ms=Omega_316;
%             c_ms=c_ms./k_ms'; 
%             Rec_ms=real(c_ms);
%             Imc_ms=imag(c_ms); 
% %             Rec_ms=Omega/Rek_ms;
% %             Imk_ms=Omega/Imk_ms; 
% clear data; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(showplot==1)
  xlab='k_0 L/\pi, dimensionless frequency';
  xlimit1=[min(klpi) 5]; %%limit is 0.05~5 
  xtick1=[0.05; 0.25; 0.5; 2.5; 5];
  xlimit2=[min(klpi) 2]; %%limit is 0.05~5 
  xtick2=[0.05; 0.25; 0.5; 2];
  xlimit3=[min(klpi) 3]; %%limit is 0.05~5 
  xtick3=[0.05; 0.25; 0.5; 2; 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cfig=figure(); %%%Re[k]
%       semilogx(klpi,Rek,'k','LineWidth',2); 
%       set(gca,'XLim',xlimit1,'XTick',xtick1,'XLabel',xlab); 
%       set(gca,'YLabel','wavevector[1/m]');
%       h=legend('Re[k]'); set(h,'Location','northeast');
%       savefigname=strcat('Rek_tabk_',filebase);
%       if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
%       if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
% 
%    cfig=figure(); %%%Im[k]
%       semilogx(klpi,Imk,'k','LineWidth',2); 
%       set(gca,'XLim',xlimit1,'XTick',xtick1,'XLabel',xlab); 
%       set(gca,'YLabel','wavevector[1/m]');
%       h=legend('Im[k]'); set(h,'Location','northeast');
%       savefigname=strcat('Imk_tabk_',filebase);
%       if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
%       if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x range from 0.05 ~ 3.0 
cfig = figure();
        semilogx(klpi,Rec,'ko','LineWidth',2);
        %hold on; semilogx(klpi_ms,real(c1),'b',klpi_ms,real(c2),'g--','LineWidth',2); 
        hold on; semilogx(klpi_ms,real(c1),'b','LineWidth',2); 
        set(gca,'XLim',xlimit3,'XTick',xtick3);
        xlabel(xlab);
        ylabel('sound speed [m/s]');
        h=legend('Re[c]'); set(h,'Location','northeast');
        savefigname=strcat('Rec1_5_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

    cfig = figure();
        semilogx(klpi,Imc,'ko','LineWidth',2);
        %hold on; semilogx(klpi_ms,imag(c1),'b',klpi_ms,imag(c2),'g--','LineWidth',2); 
        hold on; semilogx(klpi_ms,imag(c1),'b','LineWidth',2); 
        set(gca,'XLim',xlimit3,'XTick',xtick3);
        xlabel(xlab);
        ylabel('sound speed [m/s]');
        h=legend('Im[c]'); set(h,'Location','northeast');
        savefigname=strcat('Imc1_5_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x range from 0.05 ~ 2.0 
cfig = figure();
        semilogx(klpi,Rec,'ko','LineWidth',2);
        hold on; semilogx(klpi_ms,real(c1),'b',klpi_ms,real(c2),'g--','LineWidth',2); 
        set(gca,'XLim',xlimit2,'XTick',xtick2);
        xlabel(xlab);
        ylabel('sound speed [m/s]');
        h=legend('Re[c]'); set(h,'Location','northeast');
        savefigname=strcat('Rec1_2_',filebase);
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

    cfig = figure();
        semilogx(klpi,Imc,'ko','LineWidth',2);
        hold on; semilogx(klpi_ms,imag(c1),'b',klpi_ms,imag(c2),'g--','LineWidth',2); 
        set(gca,'XLim',xlimit2,'XTick',xtick2);
        xlabel(xlab);
        ylabel('sound speed [m/s]');
        h=legend('Im[c]'); set(h,'Location','northeast');
        savefigname=strcat('Imc1_2_',filebase);
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