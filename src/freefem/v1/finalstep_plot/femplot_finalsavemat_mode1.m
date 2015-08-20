%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Matlab code to analyze and plot results from FreeFem    %%%%%%%
%%%%%%% Yoonkyung Eunnie Lee 
%%%%%%% rho, chi, lambda_effective: save matlab file 
%%%%%%% Last updated on 2015-06-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set input directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basedir = '/home/eunnie12/Work/AcousticEigen_Run';
cd(basedir); 
cd('freefem/Mode1_plot/');
filebase = 'NewtCyl_Mode1_phi90_err005_finalstep';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=293.; 			   % [K]
c0=331.+0.6*(T0-273.); % [m/s]
L=1e-5;                % [m]

%---air properties
rho00=1.2;              % air density [kg/m^3]	
nuh=2.15e-5;            % kappa/(rho0*Cp)
nuv=1.5e-5;             % kinematic viscosity_v (not used) [m^2/s]
Cp=1005.;               % Heat Capacity [Joule/kg K]
kappa= 0.025929;        % Thermal conductivity, [W/mK] 
eta=1.8e-05;            % eta=nuv*rho00; [kg/ms]
Pr=0.697674418604651;	% Prandtle # = Cp eta/k , viscous/thermal diffusion
Beta0=1./T0;            % Coeff. of thermal expansion [1/K]  
gam=1.4;                % Cp/Cv 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load FEM data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(filebase,'_rhochi.txt'); data=load(filename);
klpi=data(:,1); %%% Omega is in the order of 2*pi*10^(-6)*/lambda.
xlength=length(klpi);
%klpi_100=logspace(log10(0.05),log10(5),100);
Omega = klpi*c0*pi/L;
rho=data(:,2)+data(:,3).*1i; 
chi=data(:,4)+data(:,5).*1i; 
Kbulk = 1./chi; 
filename = strcat(filebase,'_kc.txt'); data=load(filename);
k=data(:,2)+data(:,3).*1i; 
c=data(:,4)+data(:,5).*1i; 
filename = strcat(filebase,'_F.txt'); data=load(filename); 
F=data(:,2)+data(:,3).*1i; 
clear data; 
deltav= sqrt(2*eta./rho./Omega);
deltat= sqrt(2*kappa/Cp./rho./Omega);

matfilename=strcat(filebase,'.mat');
save(matfilename,'klpi','Omega','rho','chi','Kbulk','k','c','F','deltav','deltat'); 