%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
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
load('Navid_Phi90_DATA_PTS1000.mat');
data=dlmread('tabk_InitialGuess.txt'); 
tabk2=data(:,2)+i.*data(:,3);
tabk2_to_ms1=tabk2*5/6; %% this was correct 
tabk2_to_ms2=tabk2*5/4; 

klpi_100 = data(:,1); 
Omega_100 = data(:,1)*c0*pi/L; 
klpi = Vec_kalpi;
Omega = klpi*c0*pi/L;
c2_2 = Omega./K2./1e5; 
c2_tabk = Omega_100./tabk2; 
c2_tabk_1 = Omega_100./tabk2_to_ms1; 
c2_tabk_2 = Omega_100./tabk2_to_ms2; 

c1_2 = Omega./K1./1e5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(showplot==1)
  xlab='k_0 L/\pi, dimensionless frequency';
  xlimit1=[min(klpi) max(klpi)]; %%limit is 0.05~5 
  xtick1=[0.05; 0.25; 0.5; 2.5; 5];
  xlimit2=[min(klpi) 2]; %%limit is 0.05~5 
  xtick2=[0.05; 0.25; 0.5; 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   cfig=figure();
      semilogx(klpi,imag(c2_2),'r',klpi,imag(c2),'b--',klpi_100,imag(c2_tabk),'ko',klpi_100,imag(c2_tabk_1),'bo',klpi_100,imag(c2_tabk_2),'ro','LineWidth',2); 
      set(gca,'XLim',xlimit2,'XTick',xtick2); 
%      set(gca,'YLabel','sound speed [m/s]');
        %h=legend('Im[c2]'); set(h,'Location','northeast');
        savefigname='Imc2_MS';
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   cfig=figure(); 
      semilogx(klpi,real(c2_2),'r',klpi,real(c2),'b--',klpi_100,real(c2_tabk),'ko',klpi_100,real(c2_tabk_1),'bo',klpi_100,real(c2_tabk_2),'ro','LineWidth',2) 
      set(gca,'XLim',xlimit2,'XTick',xtick2); 
%      set(gca,'YLabel','sound speed [m/s]');
       % h=legend('Re[c2]'); set(h,'Location','northeast');
        savefigname='Rec2_MS';
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cfig=figure(); 
%       set(gca,'XLim',xlimit2,'XTick',xtick2); 
%  %     set(gca,'YLabel','sound speed [m/s]');
%         h=legend('Im[c]'); set(h,'Location','northeast');
%         savefigname='Imc1_MS';
%         if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
%         if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cfig=figure(); 
%       semilogx(klpi,real(c1),'b','LineWidth',2);
%       set(gca,'XLim',xlimit2,'XTick',xtick2); 
%   %    set(gca,'YLabel','sound speed [m/s]');
%         h=legend('Re[c]'); set(h,'Location','northeast');
%         savefigname='Rec1_MS';
%         if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
%         if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cfig=figure(); 
%       semilogx(klpi,imag(c3),'b','LineWidth',2);
%       set(gca,'XLim',xlimit2,'XTick',xtick2); 
%  %     set(gca,'YLabel','sound speed [m/s]');
%         h=legend('Im[c3]'); set(h,'Location','northeast');
%         savefigname='Imc3_MS';
%         if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
%         if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cfig=figure(); 
%       semilogx(klpi,real(c3),'b','LineWidth',2);
%       set(gca,'XLim',xlimit2,'XTick',xtick2); 
%   %    set(gca,'YLabel','sound speed [m/s]');
%         h=legend('Re[c3]'); set(h,'Location','northeast');
%         savefigname='Rec3_MS';
%         if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
%         if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 



end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Save as separate TXT File if asked 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(savetxt==1)
   outname= strcat('tabk_to_ms_c.txt'); 
   dlmwrite(outname,[klpi_100,real(c2_tabk_1),imag(c2_tabk_1)],'delimiter','\t','precision',prec);
   outname= strcat('tabk_to_ms_k.txt'); 
   dlmwrite(outname,[klpi_100,real(tabk2_to_ms1),imag(tabk2_to_ms1)],'delimiter','\t','precision',prec);
end %%if(savetxt) 


cd (basedir); 