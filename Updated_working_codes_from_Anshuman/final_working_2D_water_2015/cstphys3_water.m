%cstphys



%parametres air

%parametres independants ajustables
T0=20+273; %(K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho0=9.9778E+2; %#(kg/m3) Nick's experiment
%rho0=1.2; %(kg/m3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%eta=1.836924704747684e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kappa=1e-13; %0.60475; %1e-13;
kappa=0.60475; %1e-13;
%kappa=0.001; %(W/m.K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0.5, 0.3, 0.1, 0.05, O.01
%kappa=kappa*0.01;
%kappa=1e-13;
beta0=3.4112E-3; %(1/K)
%Pr=0.713;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamm=1.01;
%gamm=1.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parametres dependants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu=1*9.7937E-7; %#(m2/s)  1.5244e-5
eta=nu*rho0;

zeta=3.4e-3;%0.6.*eta;
%nu=eta/rho0; %(m2/s)  1.5244e-5
%eta=1.836924704747684e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cp=0.713*kappa/rho0/nu;  % (J/kg.K)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp=4.0764E+3;
%cp=1.0169e+03;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pr=cp*eta./kappa;
cv= cp/gamm;       %  (J/kg.K)     
p0=rho0*cp*(gamm-1.)/(gamm*beta0^2*T0);   % (Pa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0=1485.0;
%c0=sqrt(chi0inv/rho0); %  (m/s)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
chi0=1.0/rho0/c0^2; %# (1/Pa)
chi0inv=1.0/chi0;
%chi0inv=gamm.*p0;    %(Pa)          
%chi0=1./gamm/p0;    % (1/Pa)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ka=chi0inv;



Zc0=rho0*c0;         %               
cp0=sqrt(p0/rho0);   %
nup=kappa/rho0/cv; 
nupp=(zeta+eta/3)/rho0;


%parametres geometriques

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lx=9.2e-3;
%Ly=9.14e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L=1.e-2;
ht=0.2*L;  % height of main tube (m)
%ht=0.19*L;
Lx=L;
Ly=L;
Sigma=ht;
hn=0.15*L;     % height of neck (m)
%hn=0.03*L;
l=hn;
wn=0.00300*L;    % width of neck (m)
%wn=0.018*L;
sigma=wn;
wc=L-hn;       % width of cavity (m)
hc=L-2*hn-ht;


ll=l;



phi= (ht.*L + hn.*wn + hc.*wc)./(L.^2);

