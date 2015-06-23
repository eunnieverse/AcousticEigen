%cstphys



%parametres air

%parametres independants ajustables
T0=20+273; %(K)
rho0=1.2; %(kg/m3)
eta=1.836924704747684e-5;
%eta=1.836924704747684e-4;
zeta=0.6.*eta;
kappa=0.0262; %(W/m.K)
%0.5, 0.3, 0.1, 0.05, O.01
%kappa=kappa*0.01;
%kappa=1e-13;
beta0=3.43e-3; %(1/K)
%Pr=0.713;

gamm=1.4;

%parametres dependants
nu=eta/rho0; %(m2/s)  1.5244e-5
%cp=0.713*kappa/rho0/nu;  % (J/kg.K)  
cp=1.0169e+03;
Pr=cp*eta./kappa;
cv= cp/gamm;       %  (J/kg.K)     
p0=rho0*cp*(gamm-1.)/(gamm*beta0^2*T0);   % (Pa)    
chi0inv=gamm.*p0;    %(Pa)          
chi0=1./gamm/p0;    % (1/Pa)   
Ka=chi0inv;
c0=sqrt(chi0inv/rho0); %  (m/s)      
Zc0=rho0*c0;         %               
cp0=sqrt(p0/rho0);   %
nup=kappa/rho0/cv; 
nupp=(zeta+eta/3)/rho0;


%parametres geometriques

L=1.e-2;
ht=0.2*L;      % height of main tube (m)
%ht=0.19*L;
Sigma=ht;
hn=0.15*L;     % height of neck (m)
%hn=0.03*L;
l=hn;
wn=0.015*L;    % width of neck (m)
%wn=0.018*L;
sigma=wn;
wc=L-hn;       % width of cavity (m)
hc=L-2*hn-ht;


ll=l;



phi= (ht.*L + hn.*wn + hc.*wc)./(L.^2);

