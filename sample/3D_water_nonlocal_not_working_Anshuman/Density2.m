% Compute Nonlocal density - Array of Helmholtz resonators
% Using Zwikker-Konten model for slits
% Taking into account Zwikker-Konten waves in the cavity-- Horozontally



function [rho, chi]=Density2(omega,k)
 f0=10i;%-1i*k;

%constantes
%cstphys3_water;
cstphys3

rhot = rho0./Fr(nu,omega, dt/2, ht/2);         % effective density of the main tube (slit)          
rhon=  rho0./Fc(nu,omega);         % effective density of the neck (slit)
rhoc= rho0./Fr(nu,omega, dc/2, hc/2);        % effective density of the cavity (slit)


Kt= 1./chi0; %(chi0* (gamm - (gamm-1) * Fr(nup/gamm,omega, dt/2, ht/2) )   ).^(-1) ; % effective modulus of the tube (slit)
Kn= 1./chi0;% (chi0* (gamm - (gamm-1) * Fc(nup/gamm,omega) )   ).^(-1) ;               % effective modulus of the neck (slit)
Kc= 1./chi0;% (chi0* (gamm - (gamm-1) * Fr(nup/gamm,omega, dc/2, hc/2) )   ).^(-1) ;             % effective modulus of the cavity (slit)

chit=1./Kt;
chin=1./Kn;
chic=1./Kc;

ct= sqrt(Kt./rhot);             % phase velocity in the tube
cn= sqrt(Kn./rhon);             % phase velocity in the neck
cc= sqrt(Kc./rhoc);             %   phase velocity in the cavity  

kt= omega./ct;              % wave number of the tube (slit)
kn= omega./cn;              % wave number of the neck (slit)
kc= omega./cc;              % wave number of the cavity (slit)

Zct = rhot.*ct./(ht*dt);             % characteristic impedance of the main tube
Zcn= rhon.*cn./(pi*wn^2 /4);              % characteristic impedance of the neck
Zcc= rhoc.*cc./(hc*dc);            % characteristic impedance of the cavity

Yct= 1./Zct;                    % characteristic admittance of the main tube           
Ycn= 1./Zcn;                    % characteristic admittance of the neck
Ycc= 1./Zcc;                    % characteristic admittance of the cavity

% Amplitude of forcing term
%f0=10*1i;

% Constants of nonhomogeneity of equations:

% For the tube (density)
%checked
Btr= 1i.*k./(omega.^2.*rhot.*chit-k.^2);
%Btr=1i.*k./(omega.^2.*rhot./Kt-k.^2);
%checked
Ctr= (1i.*omega.*Sigma.*chit)./(omega.^2.*rhot.*chit-k.^2);
%Ctr=(1i.*omega.*Sigma./Kt)./(omega.^2.*rhot./Kt-k.^2);

% For the tube (compressibility)
Btc= (omega.^2.*rhot./(k.^2 - omega.^2.*rhot.*chit)).*(chit-gamm.*chi0);
%Btc=(1.-Kt./p0).*(-omega.^2.*rhot/Kt)./(omega.^2.*rhot/Kt-k.^2);
Ctc= (omega.*k.*Sigma./(k.^2 - omega.^2.*rhot.*chit)).*(chit-gamm.*chi0);
%Ctc=1i.*k.*(1.-Kt./p0).*(1i.*omega.*Sigma./Kt)./(omega.^2.*rhot./Kt-k.^2);

% For the neck (compressibility) -- Eqs are homog for "density" in the neck 
%Bnc= (gamm.*chi0./chin - 1).* 2./(k.*Sigma).*sin(k.*Sigma./2);

% new expression corrected for 3D
Bnc= (gamm.*chi0./chin - 1).* ( integral(@(x) 8*x.*besselj(0,k.*x)./wn^2, 0, wn/2));

%Bnc=-(1.-Kn./p0).*(2./k./sigma).*sin(k.*sigma./2);
%Ccn= 0;

% For the cavity (density)
Bcr= 1i.*k./(omega.^2.*rhoc.*chic-k.^2);%checked
%Bcr=1i.*k./(omega.^2.*rhoc./Kc-k.^2);
%Ccr= (1i.*omega.*(L-Sigma-2.*l).*chic)./(omega.^2.*rhoc.*chic-k.^2);
Ccr=(1i.*omega.*(hc)./Kc)./(omega.^2.*rhoc./Kc-k.^2);%checked

% For the cavity (compressibility)
Bcc= (omega.^2.*rhoc./(k.^2 - omega.^2.*rhoc.*chic)).*(chic-gamm.*chi0);
%Bcc=(1.-Kc./p0).*(-omega.^2.*rhoc/Kc)./(omega.^2.*rhoc/Kc-k.^2);
Ccc= (omega.*k.*(hc)./(k.^2 - omega.^2.*rhoc.*chic)).*(chic-gamm.*chi0);
%Ccc=1i.*k.*(1.-Kc./p0).*(1i.*omega.*(L-Sigma-2*l)./Kc)./(omega.^2.*rhoc./Kc-k.^2);


% Constructing the coefficient Matrix of 
%continuity relations on amplitudes



%M(:,:)=zeros(10,10);
%Crho(:,1)=zeros(10,1);
%Cchi(:,1)=zeros(10,1);


for nr=1:10
%    A(nr,1)=0.+0.*1i;
    Crho(nr,1)=0.+0.*1i;
    Cchi(nr,1)=0.+0.*1i;
    for nc=1:10
        M(nr,nc)=0.+0.*1i;
    end
end
%checked
B1= [-exp(-1i.*kt.*Lx./2).*exp(1i.*k.*Lx), -exp(1i.*kt.*Lx./2).*exp(1i.*k.*Lx), exp(1i.*kt.*Lx./2)...
    exp(-1i.*kt.*Lx./2); 
-exp(-1i.*kt.*Lx./2).*exp(1i.*k.*Lx), exp(1i.*kt.*Lx./2).*exp(1i.*k.*Lx), exp(1i.*kt.*Lx./2),...
-exp(-1i.*kt.*Lx./2);
-1, -1, 1, 1;
-1, -1, 0, 0;
Yct, -Yct, -Yct, Yct;
];

B2= zeros(3,6);

B3= zeros(5,4);
%checked
B4= [exp(-1i.*kn.*hn./2), exp(1i.*kn.*hn./2) 0, 0, 0, 0;
-Ycn.*exp(-1i.*kn.*hn./2), Ycn.*exp(1i.*kn.*hn./2), 0, 0, 0, 0;
exp(1i.*kn.*hn./2), exp(-1i.*kn.*hn./2) -1, -1, 0, 0;
Ycn.*exp(1i.*kn.*hn./2), -Ycn.*exp(-1i.*kn.*hn./2), Ycc, -Ycc, -Ycc, Ycc];
%checked
B5= [0, 0, 1, 1, -1, -1;
0, 0, Ycc.*exp((-1./2).*1i.*kc.*(wc)), -Ycc.*exp((1./2).*1i.*kc.*(wc)), 0, 0;
0, 0, 0, 0, Ycc.*exp((1./2).*1i.*kc.*(wc)), -Ycc.*exp((-1./2).*1i.*kc.*(wc))];


%checked
M= [[B1; B3] [B2; B4; B5]]; 



% Constructing the coef matrix AA by its blocks
%checked
Crho=[0; 0; 0; Btr; 0; Bcr; 0; 0; -Ccr.*exp((-1./2).*1i.*k.*(wc)); -Ccr.*exp((1./2).*1i.*k.*(wc))];      % vector of constants: AA * X = CC

Cchi=[0; 0; 0; Btc-Bnc; 0; Bcc-Bnc; 0; 0; -Ccc.*exp((-1./2).*1i.*k.*(wc)); -Ccc.*exp((1./2).*1i.*k.*(wc))];
%Solve the system of linear equations to get amplitudes in the matrix Ar
% and Ac

 
Ar= M\Crho;    % Amplitude for density
Ac= M\Cchi;



% copied and modified from average_vc
average_vr= Ctr.* ( exp(1i.*k.*Lx./2)-exp(-1i.*k.*Lx./2) ) ./1i./k +...
    Ccr.*( exp(1i.*k.*(wc)./2)-exp(-1i.*k.*(wc)./2) ) ./1i./k +...
    Yct.*Ar(1).*(1.-exp(-1i.*kt.*Lx./2))./1i./kt -Yct.*Ar(2).*(1.-exp(1i.*kt.*Lx./2))./(-1i)./kt +...
    Ycc.*Ar(7).*(1.-exp(-1i.*kc.*(wc)./2))./1i./kc -Ycc.*Ar(8).*(1.-exp(1i.*kc.*(wc)./2))./(-1i)./kc +...
    Yct.*Ar(3).*(exp(1i.*kt.*Lx./2)-1.)./1i./kt -Yct.*Ar(4).*(exp(-1i.*kt.*Lx./2)-1.)./(-1i)./kt +...
    Ycc.*Ar(9).*(exp(1i.*kc.*(wc)./2)-1.)./1i./kc -Ycc.*Ar(10).*(exp(-1i.*kc.*(wc)./2)-1.)./(-1i)./kc; 

average_vr=f0*average_vr./(Lx*Ly*Lz);





%copied from average_pvc
average_pvr=Btr.*Ctr.*(1.-exp(-1i.*k.*Lx))./(2.*1i.*k) +...
    Ar(1).*Ctr.*(1.-exp(-1i.*(k+kt).*Lx./2))./(1i.*(k+kt)) +...
    Ar(2).*Ctr.*(1.-exp(-1i.*(k-kt).*Lx./2))./(1i.*(k-kt)) +...
    Btr.*Yct.*Ar(1).*(1.-exp(-1i.*(k+kt).*Lx./2))./(1i.*(k+kt)) +...
    Yct.*Ar(1).^2.*(1.-exp(-1i.*kt.*Lx))./(2.*1i.*kt) -...
    Btr.*Yct.*Ar(2).*(1.-exp(-1i*(k-kt).*Lx./2))./(1i.*(k-kt)) -...
    Yct.*Ar(2).^2.*(1.-exp(1i.*kt.*Lx))./(-2.*1i.*kt);

average_pvr=average_pvr +Bcr.*Ccr.*(1.-exp(-1i.*k.*(wc)))./(2.*1i.*k)...
    +Ar(7).*Ccr.*(1.-exp(-1i.*(k+kc).*(wc)./2))./(1i.*(k+kc)) +...
    Ar(8).*Ccr.*(1.-exp(-1i.*(k-kc).*(wc)./2))./(1i.*(k-kc)) +...
    Bcr.*Ycc.*Ar(7).*(1.-exp(-1i.*(k+kc).*(wc)./2))./(1i.*(k+kc)) +...
    Ycc.*Ar(7).^2.*(1.-exp(-1i.*kc.*(wc)))./(2.*1i.*kc) -...
    Bcr.*Ycc.*Ar(8).*(1.-exp(-1i*(k-kc).*(wc)./2))./(1i.*(k-kc)) -...
    Ycc.*Ar(8).^2.*(1.-exp(1i.*kc.*(wc)))./(-2.*1i.*kc);








average_pvr=average_pvr +Btr.*Ctr.*(exp(1i.*k.*Lx)-1.)./(2.*1i.*k) +...
    Btr.*Yct.*Ar(3).*(exp(1i.*(k+kt).*Lx./2)-1.)./(1i.*(k+kt)) -...
    Btr.*Yct.*Ar(4).*(exp(1i.*(k-kt).*Lx./2)-1.)./(1i.*(k-kt)) +...
    Ac(3).*Ctr.*(exp(1i.*(k+kt).*Lx./2)-1.)./(1i.*(k+kt)) +...
    Yct.*Ar(3).^2.*(exp(1i.*kt.*Lx)-1.)./(2.*1i.*kt) + ...
    Ar(4).*Ctr.*(exp(1i.*(k-kt).*Lx./2)-1.)./(1i.*(k-kt)) -...
    Yct.*Ar(4).^2.*(exp(-1i.*kt.*Lx)-1.)./(-2.*1i.*kt);

average_pvr= average_pvr+Bcr.*Ccr.*(exp(1i.*k.*(wc))-1.)./(2.*1i.*k) +...
    Bcr.*Ycc.*Ar(9).*(exp(1i.*(k+kc).*(wc)./2)-1.)./(1i.*(k+kc)) -...
    Bcr.*Ycc.*Ar(10).*(exp(1i.*(k-kc).*(wc)./2)-1.)./(1i.*(k-kc)) +...
    Ar(9).*Ccr.*(exp(1i.*(k+kc).*(wc)./2)-1.)./(1i.*(k+kc)) +...
    Ycc.*Ar(9).^2.*(exp(1i.*kc.*(wc))-1.)./(2.*1i.*kc) + ...
    Ar(10).*Ccr.*(exp(1i.*(k-kc).*(wc)./2)-1.)./(1i.*(k-kc)) -...
    Ycc.*Ar(10).^2.*(exp(-1i.*kc.*(wc))-1.)./(-2.*1i.*kc);

average_pvr=f0*f0*average_pvr./(Lx*Ly*Lz);




average_vc= Ctc.* ( exp(1i.*k.*Lx./2)-exp(-1i.*k.*Lx./2) ) ./1i./k +...
    Ccc.*( exp(1i.*k.*(wc)./2)-exp(-1i.*k.*(wc)./2) ) ./1i./k +...
    Yct.*Ac(1).*(1.-exp(-1i.*kt.*Lx./2))./1i./kt -Yct.*Ac(2).*(1.-exp(1i.*kt.*Lx./2))./(-1i)./kt +...
    Ycc.*Ac(7).*(1.-exp(-1i.*kc.*(wc)./2))./1i./kc -Ycc.*Ac(8).*(1.-exp(1i.*kc.*(wc)./2))./(-1i)./kc +...
    Yct.*Ac(3).*(exp(1i.*kt.*Lx./2)-1.)./1i./kt -Yct.*Ac(4).*(exp(-1i.*kt.*Lx./2)-1.)./(-1i)./kt +...
    Ycc.*Ac(9).*(exp(1i.*kc.*(wc)./2)-1.)./1i./kc -Ycc.*Ac(10).*(exp(-1i.*kc.*(wc)./2)-1.)./(-1i)./kc;  





average_vc=average_vc./(Lx*Ly*Lz);







average_pvc=Btc.*Ctc.*(1.-exp(-1i.*k.*Lx))./(2.*1i.*k) +...
    Ac(1).*Ctc.*(1.-exp(-1i.*(k+kt).*Lx./2))./(1i.*(k+kt)) +...
    Ac(2).*Ctc.*(1.-exp(-1i.*(k-kt).*Lx./2))./(1i.*(k-kt)) +...
    Btc.*Yct.*Ac(1).*(1.-exp(-1i.*(k+kt).*Lx./2))./(1i.*(k+kt)) +...
    Yct.*Ac(1).^2.*(1.-exp(-1i.*kt.*Lx))./(2.*1i.*kt) -...
    Btc.*Yct.*Ac(2).*(1.-exp(-1i*(k-kt).*Lx./2))./(1i.*(k-kt)) -...
    Yct.*Ac(2).^2.*(1.-exp(1i.*kt.*Lx))./(-2.*1i.*kt);
average_pvc=average_pvc +Bcc.*Ccc.*(1.-exp(-1i.*k.*(wc)))./(2.*1i.*k)...
    +Ac(7).*Ccc.*(1.-exp(-1i.*(k+kc).*(wc)./2))./(1i.*(k+kc)) +...
    Ac(8).*Ccc.*(1.-exp(-1i.*(k-kc).*(wc)./2))./(1i.*(k-kc)) +...
    Bcc.*Ycc.*Ac(7).*(1.-exp(-1i.*(k+kc).*(wc)./2))./(1i.*(k+kc)) +...
    Ycc.*Ac(7).^2.*(1.-exp(-1i.*kc.*(wc)))./(2.*1i.*kc) -...
    Bcc.*Ycc.*Ac(8).*(1.-exp(-1i*(k-kc).*(wc)./2))./(1i.*(k-kc)) -...
    Ycc.*Ac(8).^2.*(1.-exp(1i.*kc.*(wc)))./(-2.*1i.*kc);

average_pvc=average_pvc +Btc.*Ctc.*(exp(1i.*k.*Lx)-1.)./(2.*1i.*k) +Btc.*Yct.*Ac(3).*(exp(1i.*(k+kt).*Lx./2)-1.)./(1i.*(k+kt)) -...
    Btc.*Yct.*Ac(4).*(exp(1i.*(k-kt).*Lx./2)-1.)./(1i.*(k-kt)) +Ac(3).*Ctc.*(exp(1i.*(k+kt).*Lx./2)-1.)./(1i.*(k+kt))...
    +Yct.*Ac(3).^2.*(exp(1i.*kt.*Lx)-1.)./(2.*1i.*kt) + Ac(4).*Ctc.*(exp(1i.*(k-kt).*Lx./2)-1.)./(1i.*(k-kt)) -...
    Yct.*Ac(4).^2.*(exp(-1i.*kt.*Lx)-1.)./(-2.*1i.*kt);
average_pvc= average_pvc+Bcc.*Ccc.*(exp(1i.*k.*(wc))-1.)./(2.*1i.*k) +Bcc.*Ycc.*Ac(9).*(exp(1i.*(k+kc).*(wc)./2)-1.)./(1i.*(k+kc)) ...
    -Bcc.*Ycc.*Ac(10).*(exp(1i.*(k-kc).*(wc)./2)-1.)./(1i.*(k-kc)) +Ac(9).*Ccc.*(exp(1i.*(k+kc).*(wc)./2)-1.)./(1i.*(k+kc))...
    +Ycc.*Ac(9).^2.*(exp(1i.*kc.*(wc))-1.)./(2.*1i.*kc) + Ac(10).*Ccc.*(exp(1i.*(k-kc).*(wc)./2)-1.)./(1i.*(k-kc)) ...
    -Ycc.*Ac(10).^2.*(exp(-1i.*kc.*(wc))-1.)./(-2.*1i.*kc);
average_pvc=average_pvc./(Lx*Ly*Lz);


average_bc=2.*1i.*Ctc.*sin(kt.*Lx./2) -Yct.*Ac(1,1).*exp(-1i.*kt.*Lx./2)...
    +Yct.*Ac(2,1).*exp(1i.*kt.*Lx./2) +Yct.*Ac(3,1).*exp(1i.*kt.*Lx./2) -Yct.*Ac(4,1).*exp(-1i.*kt.*Lx./2);
average_bc=average_bc./(1i.*omega.*Lx*Ly*Lz);


% Compute H for density and compressibility
Hr= average_pvr./average_vr;
Hc= average_pvc./average_vc;


% Compute effective nonlocal density and compressibility

%rho= phi.*(f0-1i.*k.*Hr)./(-1i.*omega.*average_vr);
%rho=phi.*(1.-1i.*k.*Hr)./(-1i.*omega.*average_vr);
rho=phi.*(f0-1i*k.*Hr)./(-1i.*omega.*average_vr);

chiinv=chi0inv.*(1.+Hc)./(average_bc.*(chi0inv./phi)+gamm);
chi=1./chiinv;

%rhoaf=porosite.*(1.-1i.*k.*Hrho)./(-1i.*omega.*vmoyrho);
































