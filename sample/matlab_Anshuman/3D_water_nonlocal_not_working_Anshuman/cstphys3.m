%cstphys_Water



%parametres water

%parametres independants ajustables
T0=20+273; %(K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho0=1*1e3; %#(kg/m3) Nick's experiment
kappa=1e-13;%0.60475; %1e-13;
%beta0=2e-4;%3.4112E-3; %(1/K)
gamm=1.01;

%parametres dependants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu=1*9.7937E-7; %#(m2/s)  1.5244e-5
%eta=nu*rho0;
%zeta=3.4e-3;%0.6.*eta;
%cp=4.0764E+3;
%Pr=cp*eta./kappa;
%cv= cp/gamm;       %  (J/kg.K)     
%p0=rho0*cp*(gamm-1.)/(gamm*beta0^2*T0);   % (Pa)

c0=1485.0;
chi0=1.0/rho0/c0^2; %# (1/Pa)
chi0inv=1.0/chi0;
%Ka=chi0inv;
%Zc0=rho0*c0;         %               
%cp0=sqrt(p0/rho0);   %
%nup=kappa/rho0/cv; 
%nupp=(zeta+eta/3)/rho0;






%numerical values from Nick's paper
factor=1;

wc=5.0e-3 / factor;
hc=3.14e-3 / factor;
dc=4e-3 / factor;
hn=1*1e-3 / factor;
wn=1*1e-3 / factor;
xn=4.2e-3 / factor; %xn/2 is the side
ht=1*4e-3 / factor;
zn=1e-3 /factor; % half length
% derived quantities
yn=hn; % full length
Lx=wc + xn/2 + xn/2;
Ly=hc + yn + hn + ht;
Lz=dc + zn/2 + zn/2;
dt=ht;
rw=wn / 2;
l=xn;
sigma= wn;
Sigma= ht;
phi= (ht.*Lx*dt + hn.*pi*wn^2 / 4 + hc.*wc*dc)./(Lx*Ly*Lz);


% 
% Lx=9.2e-3;
% ht=4e-3;      % height of main tube (m)
% Sigma=ht;
% hn=1e-3; % height of neck (m)
% xn=4.2e-3;
% l=xn;
% wn=1*1e-3;    % width of neck (m)
% sigma=wn;
% wc=5e-3;%5e-3;       % width of cavity (m)
% hc=3.14e-3;%5e-3;%3.14e-3; %L-2*hn-ht;
% %ln=hn;
% rw=wn/2;
% %ll=l;
% dc=4e-3;
% dt=ht; % ht


%phi= (ht.*L + hn.*wn + hc.*wc)./(L.^2);
