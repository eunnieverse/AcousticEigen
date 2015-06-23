%cstphys



%parametres air

%parametres independants ajustables
T0=20+273; %(K)
rho0=1.2; %(kg/m3)
eta=1*1.836924704747684e-5;
%eta=1.836924704747684e-4;
zeta=0.6.*eta;
kappa=0.0257;%1e-6;%0.0257; %(W/m.K)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lx=9.2e-3;
%Ly=9.14e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% real wn=0.5*1e-3;       // width of neck (m)
% real ht=0.5*4e-3;          // height of main tube (m)
% real hn=(pi*wn*(1.0e-3))/(4.0*ht);        // height of neck (m)
% 
% real Lx=9.2e-3;
% real wc=Lx-4.2e-3;       // width of cavity (m)
% real hc=(5e-3*3.14e-3*4e-3)/(4*ht*wc);  // height of cavity (m)
% real Ly=ht + hc + 1e-3 + hn;           // length of cell /y (m)


%% updated geometrical parameters (1st Nov. 2014)
Lxp=1e-2;
wnp=0.6*0.080*Lxp; 
htp=0.2*Lxp;  
hnp=0.15*Lxp; 
xnp=hnp;
wcp=Lxp-xnp;
hcp=Lxp-2*hnp-htp;
dcp=hcp;


ynp=hnp; %chosen

%3D geometry has ended

% now we start to convert these to 2D

Lx=1*Lxp; % assumption
xn=1*xnp; %assumption
yn=ynp; % assumption

wn=0.5*wnp; %width of neck (m)
ht=0.5*htp; %height of main tube (m)
%hn=(pi*wn*hnp*Lx)/(4.0*ht*Lxp); %height of neck (m)





wc=Lx-xn;       % width of cavity (m)
%hc=(wcp*hcp*dcp*Lx)/(4*ht*wc*Lxp);  % original 2D version of 3D height of cavity (m)
hc=(hcp*wcp*dcp)/(2*htp*wc); %Navid's Jan 2015 email

 hn=(hnp*wcp*hcp*dcp*wnp/2)/((pi*wnp*wnp/4)*wc*hc);


Ly=ht + hc + yn + hn;           % length of cell Ly (m)
phi= (ht.*Lx + hn.*wn + hc.*wc)./(Lx*Ly);

%% old set of parameters
% Lx=1.e-3;
% Ly=1.e-3;
% ht=0.2*Lx;  % height of main tube (m)
% hn=0.15*Lx;     % height of neck (m)
% wn=0.01200*Lx;    % width of neck (m)
% wc=Lx-hn;       % width of cavity (m)
% hc=Ly-2*hn-ht;
% phi= (ht.*Lx + hn.*wn + hc.*wc)./(Lx*Ly);

