%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
%%%%%%% search for tabk definition section in *.edp and plot it 
%%%%%%% Last updated on 2015-06-10
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
basedir = pwd(); %% directory containing .m file 
%basedir = '/home/eunnie12/Dropbox (MIT)/labwork/Presentations_2015/11_AcousticMeta/Eunnie_Navid';
cd(basedir); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=293.; 			   % [K]
c0=331.+0.6*(T0-273.); % [m/s]
L=1e-5;                % [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read data into Rek, Imk 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirlist=dir('*.txt');     %%for the first .txt file in the folder
filename=dirlist.name; 
fnsplit = strsplit(filename,'.'); %% fn= 1x2 cell, 'abcd_efgh' 'edp' 
fnext = char(fnsplit(length(fnsplit))); %% 'edp' 
filebase = char(fnsplit(1));

fid = fopen(filename);
   formatSpec = 'tabk(%*u)=%f;     // Vec_kaLpi=%*f\n';
   data = textscan(fid,formatSpec);
   data = data{1,1};
   fclose(fid); 

%%% define x axis: nondimensional frequency klpi (kL/pi), Omega
klpi=logspace(log10(0.05),log10(5),length(data));       
Omega= klpi*c0*pi/L;
%%% define y data: 
k=data;
Rek=real(k);
Imk=imag(k);
clear data; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(showplot==1)

   cfig=figure(); %%%Re[k]
      semilogx(klpi,Rek,'k','LineWidth',2); 
      xlim([min(klpi) max(klpi)]); xlabel('k_0 L/\pi, dimensionless frequency');
      set(gca,'XTick',[0.05;0.1;0.5;1;2;5]);
      h=legend('Re[k]'); set(h,'Location','northeast');
      savefigname='Rek_tabk';
      if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
      if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

   cfig=figure(); %%%Im[k]
      semilogx(klpi,Imk,'k','LineWidth',2); 
      xlim([min(klpi) max(klpi)]); xlabel('k_0 L/\pi, dimensionless frequency');
      set(gca,'XTick',[0.05;0.1;0.5;1;2;5]);
      h=legend('Im[k]'); set(h,'Location','northeast');
      savefigname='Imk_tabk';
      if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
      if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

end %%if(showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save as separate TXT File if asked 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(savetxt==1)
   outname= strcat('tabk','_',filebase,'.txt'); 
   dlmwrite(outname,[klpi',Rek,Imk],'delimiter','\t','precision',prec);
end %%if(savetxt) 

cd (basedir); 