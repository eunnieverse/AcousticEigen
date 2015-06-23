%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
%%%%%%% Local
%%%%%%% Last updated on 2015-06-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set output plot conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extension = '.png'; %%default plot extension 
savefig = 1; 
saveeps = 0; 
savetxt = 1; 
savemat = 0; 
prec = '%e'; %% precision 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set input directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basedir = '/home/eunnie12/Work/AcousticEigen';
cd(basedir); 
cd('freefem/Local/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=293.; 			   % [K]
c0=331.+0.6*(T0-273.); % [m/s]
L=1e-5;                % [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load rhoLocal.txt
load chiLocal.txt
Omega=2.*pi.*rhoLocal(:,1);
klpi=(Omega.*L)./(c0.*pi);
cL=c0./sqrt((rhoLocal(:,4)+1i.*rhoLocal(:,5)).*(chiLocal(:,4)+1i.*chiLocal(:,5)));
Rec = real(cL);
Imc = imag(cL); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlab='kL/\pi, dimensionless frequency';
  xlimit1=[min(klpi) max(klpi)]; %%limit is 0.05~5 
  xtick1=[0.05; 0.25; 0.5; 2.5; 5]; 
cfig = figure();
        semilogx(klpi,Rec,'--r','LineWidth',2);
        set(gca,'XLim',xlimit1,'XTick',xtick1);
        xlabel(xlab);
        ylabel('phase velocity [m/s]');
        h=legend('Re[c]'); set(h,'Location','northeast');
        savefigname='Local_Rec';
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 

        cfig = figure();
        semilogx(klpi,Imc,'--r','LineWidth',2);
        set(gca,'XLim',xlimit1,'XTick',xtick1);
        xlabel(xlab);
        ylabel('phase velocity [m/s]');
        h=legend('Im[c]'); set(h,'Location','northeast');
        savefigname='Local_Imc';
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
        
 %% save local model c into txt file
if(savetxt==1)
   outname= strcat('c_local','.txt'); 
   dlmwrite(outname,[klpi,Rec,Imc],'delimiter','\t','precision',prec);
   outname= strcat('k_local','.txt'); 
   dlmwrite(outname,[klpi,Rec,Imc],'delimiter','\t','precision',prec);

end %%if(savetxt) 