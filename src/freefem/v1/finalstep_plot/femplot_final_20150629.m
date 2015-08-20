%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
%%%%%%% rho, chi, lambda_effective
%%%%%%% Last updated on 2015-06-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set output plot conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extension = '.jpg'; %%default plot extension 
savefig = 1; 
saveeps = 0; 
savetxt = 1; 
savemat = 0; 
prec = '%e'; %% precision 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set input directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basedir = '/home/eunnie12/Work/AcousticEigen_Run';
cd(basedir); 
cd('freefem/finalstep_plot/');
filebase = 'NewtCyl_phi90_err002_finalstep';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=293.; 			   % [K]
c0=331.+0.6*(T0-273.); % [m/s]
L=1e-5;                % [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load FEM data for mode 1, mode 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M1=load('NewtCyl_Mode1_phi90_err002_finalstep.mat');
M2=load('NewtCyl_Mode2_phi90_err002_finalstep.mat');

klpi=M1.klpi;

Omega=M1.Omega;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load local model data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ../Local_correct/rhoLocal.txt;
load ../Local_correct/chiLocal.txt;
load ../Local_correct/c_local.txt;
load ../Local_correct/k_local.txt;
rhoL = real(rhoLocal(:,2))+imag(rhoLocal(:,3)).*1i;
chiL = real(chiLocal(:,2))+imag(chiLocal(:,3)).*1i;
cL = real(c_local(:,2))+imag(c_local(:,3)).*1i;
kL = real(k_local(:,2))+imag(k_local(:,3)).*1i; 
Rec = real(cL); Imc = imag(cL); 
Rek = real(kL); Imk = imag(kL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlab='kL/\pi, dimensionless frequency';
  xlimit1=[min(klpi) max(klpi)]; %%limit is 0.05~5 
  xtick1=[0.05; 0.25; 0.5; 2.5; 5]; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig = figure(); %%Rho1
    semilogx(klpi,real(M1.rho),'bo',M2.klpi,real(M2.rho),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('Re(\rho)');
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylabel('effective density \rho [kg/m^3]');
    savefigname=strcat('rho1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%Rho2
    semilogx(klpi,imag(M1.rho),'bo',M2.klpi,imag(M2.rho),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('Im(\rho)');
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylabel('effective density \rho [kg/m^3]');
    savefigname=strcat('rho2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%chi1
    semilogx(klpi,real(M1.chi),'bo',M2.klpi,real(M2.chi),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('Re(\chi)');
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylabel('effective compressibility \chi [ms^2/kg]');
    savefigname=strcat('chi1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%chi2
    semilogx(klpi,imag(M1.chi),'bo',M2.klpi,imag(M2.chi),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('Im(\chi)');
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylabel('effective compressibility \chi [ms^2/kg]');
    savefigname=strcat('chi2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%K1
    semilogx(klpi,real(M1.Kbulk),'bo',M2.klpi,real(M2.Kbulk),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1);xlabel(xlab);
    title('Re(K)')
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylim([-1e05 5e05]); ylabel('bulk modulus K [kg/ms^2]');
    savefigname=strcat('Kbulk1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%K2
    semilogx(klpi,imag(M1.Kbulk),'bo',M2.klpi,imag(M2.Kbulk),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    ylim([-1e05 3e05]); ylabel('bulk modulus K [kg/ms^2]');
    title('Re(K)')
    h=legend('mode1','mode2');set(h,'Location','northeast');
    savefigname=strcat('Kbulk2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig= figure(); %% F
    semilogx(klpi,real(M1.F),'bo','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    ylabel('F=\rho\chi\omega^2-k^2 [1/m^2]');
    title('Newton F, Mode1');
    savefigname=strcat('F_mode1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig= figure(); %% F
    semilogx(M2.klpi,real(M2.F),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    ylabel('F=\rho\chi\omega^2-k^2 [1/m^2]');
    title('Newton F, Mode2');
    savefigname=strcat('F_mode2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig = figure(); %c1
    semilogx(klpi,real(M1.c),'bo',M2.klpi,real(M2.c),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    h=legend('mode1','mode2');set(h,'Location','northwest');
    title('Re(c)');
    ylabel('effective phase velocity c[m/s]');
    savefigname=strcat('c1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %c2
    semilogx(klpi,imag(M1.c),'bo',M2.klpi,imag(M2.c),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    h=legend('mode1','mode2');set(h,'Location','northeast');
    title('Im(c)');
    ylabel('effective phase velocity c[m/s]');
    savefigname=strcat('c2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig = figure(); %k1
    semilogx(klpi,real(M1.k),'bo',M2.klpi,real(M2.k),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    h=legend('mode1','mode2');set(h,'Location','northeast');
    title('Re(k)');
    ylabel('effective wavenumber k[1/m]');
    savefigname=strcat('k1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %k2
    semilogx(klpi,imag(M1.k),'bo',M2.klpi,imag(M2.k),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    h=legend('mode1','mode2');set(h,'Location','northeast');
    title('Im(k)');
    ylabel('effective wavenumber k[1/m]');
    savefigname=strcat('k2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig = figure(); %%lambda_effective
    semilogx(klpi,real(2*pi/L./M1.k),'bo','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('effective wavelength/L, Mode1');
    ylabel('effective \lambda/L');
    savefigname=strcat('lambdaeff_mode1_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%lambda_effective
    semilogx(M2.klpi,real(2*pi/L./M2.k),'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('effective wavelength/L, Mode2');
    ylabel('effective \lambda/L');
    savefigname=strcat('lambdaeff_mode2_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfig = figure(); %%deltav
    semilogx(klpi,real(M1.deltav)*1e6,'bo',M2.klpi,real(M2.deltav)*1e6,'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('viscous BL \delta_v');
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylabel('boundary layer thicnkess[\mum]');
    savefigname=strcat('deltav_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
cfig = figure(); %%deltat
    semilogx(klpi,real(M1.deltat)*1e6,'bo',M2.klpi,real(M2.deltat)*1e6,'go','LineWidth',2);
    set(gca,'XLim',xlimit1,'XTick',xtick1); xlabel(xlab);
    title('thermal BL \delta_t');
    h=legend('mode1','mode2');set(h,'Location','northeast');
    ylabel('boundary layer thicnkess[\mum]');
    savefigname=strcat('deltat_',filebase);
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
