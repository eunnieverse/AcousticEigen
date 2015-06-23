% Compute Nonlocal density - Array of Helmholtz resonators
% Using Zwikker-Konten model for slits
% Taking into account Zwikker-Konten waves in the cavity-- Horozontally



function [rho, chi]=Density2(omega,k)
%constantes
cstphys3;

st=sqrt( omega.*rho0.*((ht./2).^2)./eta );                   % define parameter
sn=sqrt( omega.*rho0.*((wn./2).^2)./eta );                   % define parameter
sc=sqrt( omega.*rho0.*((hc./2).^2)./eta );                   % define parameter


rhot = rho0./( 1-tanh(st.*sqrt(-1i))./(st.*sqrt(-1i)) );         % effective density of the main tube (slit)          
rhon=  rho0./( 1-tanh(sn.*sqrt(-1i))./(sn.*sqrt(-1i)) );         % effective density of the neck (slit)
rhoc= rho0./( 1-tanh(sc.*sqrt(-1i))./(sc.*sqrt(-1i)) );          % effective density of the cavity (slit)

Kc=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*sc.*sqrt(-1i))./(sqrt(Pr).*sc.*sqrt(-1i)) );              % effective modulus of the cavity (slit)
Kt= gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*st.*sqrt(-1i))./(sqrt(Pr).*st.*sqrt(-1i)) );               % effective modulus of the tube (slit)
Kn=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*sn.*sqrt(-1i))./(sqrt(Pr).*sn.*sqrt(-1i)) );              % effective modulus of the neck (slit)
  
chit=1./Kt;
chin=1./Kn;
chic=1./Kc;

ct= sqrt(Kt./rhot);             % phase velocity in the tube
cn= sqrt(Kn./rhon);             % phase velocity in the neck
cc= sqrt(Kc./rhoc);             %   phase velocity in the cavity  

kt= omega./ct;              % wave number of the tube (slit)
kn= omega./cn;              % wave number of the neck (slit)
kc= omega./cc;              % wave number of the cavity (slit)

Zct = rhot.*ct./ht;             % characteristic impedance of the main tube
Zcn= rhon.*cn./wn;              % characteristic impedance of the neck
Zcc= rhoc.*cc./hc;              % characteristic impedance of the cavity

Yct= 1./Zct;                    % characteristic admittance of the main tube           
Ycn= 1./Zcn;                    % characteristic admittance of the neck
Ycc= 1./Zcc;                    % characteristic admittance of the cavity

% Amplitude of forcing term
%f0=10*1i;

% Constants of nonhomogeneity of equations:

% For the tube (density)
Btr= 1i.*k./(omega.^2.*rhot.*chit-k.^2);
%Btr=1i.*k./(omega.^2.*rhot./Kt-k.^2);
Ctr= (1i.*omega.*Sigma.*chit)./(omega.^2.*rhot.*chit-k.^2);
%Ctr=(1i.*omega.*Sigma./Kt)./(omega.^2.*rhot./Kt-k.^2);

% For the tube (compressibility)
Btc= (omega.^2.*rhot./(k.^2 - omega.^2.*rhot.*chit)).*(chit-gamm.*chi0);
%Btc=(1.-Kt./p0).*(-omega.^2.*rhot/Kt)./(omega.^2.*rhot/Kt-k.^2);
Ctc= (omega.*k.*Sigma./(k.^2 - omega.^2.*rhot.*chit)).*(chit-gamm.*chi0);
%Ctc=1i.*k.*(1.-Kt./p0).*(1i.*omega.*Sigma./Kt)./(omega.^2.*rhot./Kt-k.^2);

% For the neck (compressibility) -- Eqs are homog for "density" in the neck 
Bnc= (gamm.*chi0./chin - 1).* 2./(k.*Sigma).*sin(k.*Sigma./2);
%Bnc=-(1.-Kn./p0).*(2./k./sigma).*sin(k.*sigma./2);
%Ccn= 0;

% For the cavity (density)
Bcr= 1i.*k./(omega.^2.*rhoc.*chic-k.^2);
%Bcr=1i.*k./(omega.^2.*rhoc./Kc-k.^2);
%Ccr= (1i.*omega.*(L-Sigma-2.*l).*chic)./(omega.^2.*rhoc.*chic-k.^2);
Ccr=(1i.*omega.*(L-Sigma-2*l)./Kc)./(omega.^2.*rhoc./Kc-k.^2);

% For the cavity (compressibility)
Bcc= (omega.^2.*rhoc./(k.^2 - omega.^2.*rhoc.*chic)).*(chic-gamm.*chi0);
%Bcc=(1.-Kc./p0).*(-omega.^2.*rhoc/Kc)./(omega.^2.*rhoc/Kc-k.^2);
Ccc= (omega.*k.*(L-Sigma-2.*l)./(k.^2 - omega.^2.*rhoc.*chic)).*(chic-gamm.*chi0);
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

B1= [-exp(-1i.*kt.*L./2).*exp(1i.*k.*L), -exp(1i.*kt.*L./2).*exp(1i.*k.*L), exp(1i.*kt.*L./2) exp(-1i.*kt.*L./2); 
-exp(-1i.*kt.*L./2).*exp(1i.*k.*L), exp(1i.*kt.*L./2).*exp(1i.*k.*L), exp(1i.*kt.*L./2), -exp(-1i.*kt.*L./2);
1, 1, -1, -1;
-1, -1, 0, 0;
Yct, -Yct, -Yct, Yct;
];

B2= zeros(3,6);

B3= zeros(5,4);

B4= [exp(-1i.*kn.*l./2), exp(1i.*kn.*l./2) 0, 0, 0, 0;
-Ycn.*exp(-1i.*kn.*l./2), Ycn.*exp(1i.*kn.*l./2), 0, 0, 0, 0;
exp(1i.*kn.*l./2), exp(-1i.*kn.*l./2) -1, -1, 0, 0;
Ycn.*exp(1i.*kn.*l./2), -Ycn.*exp(-1i.*kn.*l./2), Ycc, -Ycc, -Ycc, Ycc];

B5= [0, 0, 1, 1, -1, -1;
0, 0, Ycc.*exp((-1./2).*1i.*kc.*(L-l)), -Ycc.*exp((1./2).*1i.*kc.*(L-l)), 0, 0;
0, 0, 0, 0, Ycc.*exp((1./2).*1i.*kc.*(L-l)), -Ycc.*exp((-1./2).*1i.*kc.*(L-l))];



M= [[B1; B3] [B2; B4; B5]];                                              % Constructing the coef matrix AA by its blocks
Crho=[0; 0; 0; Btr; 0; Bcr; 0; 0; -Ccr.*exp((-1./2).*1i.*k.*(L-l)); -Ccr.*exp((1./2).*1i.*k.*(L-l))];      % vector of constants: AA * X = CC
Cchi=[0; 0; 0; Btc-Bnc; 0; Bcc-Bnc; 0; 0; -Ccc.*exp((-1./2).*1i.*k.*(L-l)); -Ccc.*exp((1./2).*1i.*k.*(L-l))];
%Solve the system of linear equations to get amplitudes in the matrix Ar
% and Ac

 




% M(1,1)=-exp(1i.*k*L).*exp(-1i.*kt.*L./2);
% M(1,2)=-exp(1i.*k*L).*exp(1i.*kt.*L./2);
% M(1,3)=exp(1i.*kt.*L./2);
% M(1,4)=exp(-1i.*kt.*L./2);
% 
% M(2,1)=-exp(1i.*k*L).*exp(-1i.*kt.*L./2);
% M(2,2)=exp(1i.*k*L).*exp(1i.*kt.*L./2);
% M(2,3)=exp(1i.*kt.*L./2);
% M(2,4)=-exp(-1i.*kt.*L./2);
% 
% M(3,1)=1.;
% M(3,2)=1.;
% M(3,3)=-1.;
% M(3,4)=-1.;
% 
% M(4,1)=-1.;
% M(4,2)=-1.;
% M(4,5)=exp(-1i.*kn.*l./2);
% M(4,6)=exp(1i.*kn.*l./2);
% 
% M(5,1)=Yct;
% M(5,2)=-Yct;
% M(5,3)=-Yct;
% M(5,4)=Yct;
% M(5,5)=-Ycn.*exp(-1i.*kn.*l./2);
% M(5,6)=Ycn.*exp(1i.*kn.*l./2);
% 
% M(6,5)=exp(1i.*kn.*l./2);
% M(6,6)=exp(-1i.*kn.*l./2);
% M(6,7)=-1.;
% M(6,8)=-1.;
% 
% M(7,5)=Ycn.*exp(1i.*kn.*l./2);
% M(7,6)=-Ycn.*exp(-1i.*kn.*l./2);
% M(7,7)=Ycc;
% M(7,8)=-Ycc;
% M(7,9)=-Ycc;
% M(7,10)=Ycc;
% 
% M(8,7)=1.;
% M(8,8)=1.;
% M(8,9)=-1.;
% M(8,10)=-1.;
% 
% M(9,7)=Ycc.*exp(-1i.*kc.*(L-l)./2);
% M(9,8)=-Ycc.*exp(1i.*kc.*(L-l)./2);
% 
% M(10,9)=Ycc.*exp(1i.*kc.*(L-l)./2);
% M(10,10)=-Ycc.*exp(-1i.*kc.*(L-l)./2);

% Crho(4,1)=Btr;
% Crho(6,1)=Bcr;
% Crho(9,1)=-Ccr.*exp(-1i.*k.*(L-l)./2);
% Crho(10,1)=-Ccr.*exp(1i.*k.*(L-l)./2);
% 
% Cchi(4,1)=Btc-Bnc;
% Cchi(6,1)=Bcc-Bnc;
% Crho(9,1)=-Ccc.*exp(-1i.*k.*(L-l)./2);
% Cchi(10,1)=-Ccc.*exp(1i.*k.*(L-l)./2);




Ar= M\Crho;    % Amplitude for density
Ac= M\Cchi;

%  Amplitudes are now known  

%Compute average of the velocity: <v> for density
% average_vr= f0/L^2.*((Yct.*Ar(1,1)./1i./kt).*(1 - exp(-1i.*kt.*L/2)) ...
% + (Yct.*Ar(2,1)./1i./kt).*(1 -  exp(1i.*kt.*L/2)) ...
% + (Ctr./1i./k).*(1 - exp(-1i.*k.*L/2)) ...
% + (Yct.*Ar(3,1)./1i./kt).*(exp(1i.*kt.*L/2) - 1) ...
% + (Yct.*Ar(4,1)./1i./kt).*(exp(-1i.*kt.*L/2) - 1) ...
% + (Ctr./1i./k).*(exp(1i.*k.*L/2) - 1) ...
% + (Ycc.*Ar(7,1)./1i./kc).*(1 - exp((-1/2).*1i.*kc.*(L-l))) ...
% + (Ycc.*Ar(8,1)./1i./kc).*(1 - exp((1/2).*1i.*kc.*(L-l))) ...
% + (Ccr./1i./k).*(1 - exp((-1./2).*1i.*k.*(L-l))) ...
% + (Ycc.*Ar(9,1)./1i./kc).*(exp((1/2).*1i.*kc.*(L-l)) - 1) ...
% + (Ycc.*Ar(10,1)./1i./kc).*(exp((-1./2).*1i.*kc.*(L-l)) - 1) ...
% + (Ccr./1i./k).*(exp((1./2).*1i.*k.*(L-l)) -1));

average_vr= Ctr.* ( exp(1i.*k.*L./2)-exp(-1i.*k.*L./2) ) ./1i./k +Ccr.*( exp(1i.*k.*(L-l)./2)-exp(-1i.*k.*(L-l)./2) ) ./1i./k +Yct.*Ar(1).*(1.-exp(-1i.*kt.*L./2))./1i./kt -Yct.*Ar(2).*(1.-exp(1i.*kt.*L./2))./(-1i)./kt +Ycc.*Ar(7).*(1.-exp(-1i.*kc.*(L-l)./2))./1i./kc -Ycc.*Ar(8).*(1.-exp(1i.*kc.*(L-l)./2))./(-1i)./kc + Yct.*Ar(3).*(exp(1i.*kt.*L./2)-1.)./1i./kt -Yct.*Ar(4).*(exp(-1i.*kt.*L./2)-1.)./(-1i)./kt +Ycc.*Ar(9).*(exp(1i.*kc.*(L-l)./2)-1.)./1i./kc -Ycc.*Ar(10).*(exp(-1i.*kc.*(L-l)./2)-1.)./(-1i)./kc;  
average_vr=average_vr./L^2;


%Compute <pv> for density
% average_pvr= f0/(L^2)*(((Yct*Ar(1,1)^2)/(2*1i*kt))*(1 - exp(-1i*kt*L)) ...
% + ((Ar(1,1).*Ctr)./(1i.*(kt+k))).*(1 - exp((-1./2).*1i.*(kt+k).*L)) ...
% + ((Yct.*Ar(2,1).^2)./(2.*1i.*kt)).*(1 - exp(1i.*kt.*L)) ...
% + ((Ar(2,1).*Ctr)./(1i.*(k-kt))).*(1 - exp((-1./2).*1i.*(k-kt).*L)) ...
% + ((Yct.*Btr.*Ar(1,1))./(1i.*(kt+k))).*(1 - exp((-1./2).*1i.*(kt+k).*L)) ...
% - ((Yct.*Btr.*Ar(2,1))./(1i.*(k-kt))).*(1 - exp((-1./2).*1i.*(k-kt).*L)) ...
% + ((Btr.*Ctr)./(2.*1i.*k)).*(1 - exp(-1i.*k.*L)) ...
% + ((Yct.*Ar(3,1).^2)./(2.*1i.*kt)).*(exp(1i.*kt.*L) - 1) ...
% + ((Ar(3,1).*Ctr)./(1i.*(kt+k))).*(exp((1./2).*1i.*(kt+k).*L) -1) ...
% + ((Yct.*Ar(4,1).^2)./(2.*1i.*kt)).*(exp(-1i.*kt.*L) - 1) ...
% + ((Ar(4,1).*Ctr)./(1i.*(k-kt))).*(exp((1./2).*1i.*(k-kt).*L) - 1) ...
% + ((Yct.*Btr.*Ar(3,1))./(1i.*(kt+k))).*(exp((1./2).*1i.*(kt+k).*L) - 1) ...
% - ((Yct.*Btr.*Ar(4,1))./(1i.*(k-kt))).*(exp((1./2).*1i.*(k-kt).*L) - 1) ...
% + ((Btr.*Ctr)./(2.*1i.*k)).*(exp(1i.*k.*L) - 1) ...
% + ((Ycc.*Ar(7,1).^2)./(2.*1i.*kc)).*(1 - exp(-1i.*kc.*(L-l))) ...
% + ((Ar(7,1).*Ccr)./(1i.*(kc+k))).*(1 - exp((-1./2).*1i.*(kc+k).*(L-l))) ...
% + ((Ycc.*Ar(8,1).^2)./(2.*1i.*kc)).*(1 - exp(1i.*kc.*(L-l))) ...
% + ((Ar(8,1).*Ccr)./(1i.*(k-kc))).*(1 - exp((-1./2).*1i.*(k-kc).*(L-l))) ...
% + ((Ycc.*Bcr.*Ar(7,1))./(1i.*(kc+k))).*(1 - exp((1./2).*1i.*(kc+k).*(L-l))) ...
% - ((Ycc.*Bcr.*Ar(8,1))./(1i.*(k-kc))).*(1 - exp((-1./2).*1i.*(k-kc).*(L-l))) ...
% + ((Ycc.*Ar(9,1).^2)./(2.*1i.*kc)).*(exp(1i.*kc.*(L-l)) - 1) ...
% + ((Ar(9,1).*Ccr)./(1i.*(kc+k))).*(exp((1./2).*1i.*(kc+k).*(L-l)) - 1) ...
% + ((Ycc.*Ar(10,1).^2)./(2.*1i.*kc)).*(exp(-1i.*kc.*(L-l)) - 1) ...
% + ((Ar(10,1).*Ccr)./(1i.*(k-kc))).*(exp((1./2).*1i.*(k-kc).*(L-l)) - 1) ...
% + ((Ycc.*Bcr.*Ar(9,1))./(1i.*(kc+k))).*(exp((1./2).*1i.*(kc+k).*(L-l)) - 1) ...
% - ((Ycc.*Bcr.*Ar(10,1))./(1i.*(k-kc))).*(exp((1./2).*1i.*(k-kc).*(L-l)) - 1) ...
% + ((Bcr.*Ccr)./(2.*1i.*k)).*(exp(1i.*k.*(L-l)) - 1) ...
% + ((Bcr.*Ccr)./(2.*1i.*k)).*(1 - exp(-1i.*k.*(L-l))));

average_pvr=Btr.*Ctr.*(1.-exp(-1i.*k.*L))./(2.*1i.*k) +Ar(1).*Ctr.*(1.-exp(-1i.*(k+kt).*L./2))./(1i.*(k+kt)) +Ar(2).*Ctr.*(1.-exp(-1i.*(k-kt).*L./2))./(1i.*(k-kt)) +Btr.*Yct.*Ar(1).*(1.-exp(-1i.*(k+kt).*L./2))./(1i.*(k+kt)) +Yct.*Ar(1).^2.*(1.-exp(-1i.*kt.*L))./(2.*1i.*kt) -Btr.*Yct.*Ar(2).*(1.-exp(-1i*(k-kt).*L./2))./(1i.*(k-kt)) -Yct.*Ar(2).^2.*(1.-exp(1i.*kt.*L))./(-2.*1i.*kt);
average_pvr=average_pvr +Bcr.*Ccr.*(1.-exp(-1i.*k.*(L-l)))./(2.*1i.*k) +Ar(7).*Ccr.*(1.-exp(-1i.*(k+kc).*(L-l)./2))./(1i.*(k+kc)) +Ar(8).*Ccr.*(1.-exp(-1i.*(k-kc).*(L-l)./2))./(1i.*(k-kc)) +Bcr.*Ycc.*Ar(7).*(1.-exp(-1i.*(k+kc).*(L-l)./2))./(1i.*(k+kc)) +Ycc.*Ar(7).^2.*(1.-exp(-1i.*kc.*(L-l)))./(2.*1i.*kc) -Bcr.*Ycc.*Ar(8).*(1.-exp(-1i*(k-kc).*(L-l)./2))./(1i.*(k-kc)) -Ycc.*Ar(8).^2.*(1.-exp(1i.*kc.*(L-l)))./(-2.*1i.*kc);
average_pvr=average_pvr +Btr.*Ctr.*(exp(1i.*k.*L)-1.)./(2.*1i.*k) +Btr.*Yct.*Ar(3).*(exp(1i.*(k+kt).*L./2)-1.)./(1i.*(k+kt)) -Btr.*Yct.*Ar(4).*(exp(1i.*(k-kt).*L./2)-1.)./(1i.*(k-kt)) +Ar(3).*Ctr.*(exp(1i.*(k+kt).*L./2)-1.)./(1i.*(k+kt)) +Yct.*Ar(3).^2.*(exp(1i.*kt.*L)-1.)./(2.*1i.*kt) + Ar(4).*Ctr.*(exp(1i.*(k-kt).*L./2)-1.)./(1i.*(k-kt)) -Yct.*Ar(4).^2.*(exp(-1i.*kt.*L)-1.)./(-2.*1i.*kt);
average_pvr= average_pvr+Bcr.*Ccr.*(exp(1i.*k.*(L-l))-1.)./(2.*1i.*k) +Bcr.*Ycc.*Ar(9).*(exp(1i.*(k+kc).*(L-l)./2)-1.)./(1i.*(k+kc)) -Bcr.*Ycc.*Ar(10).*(exp(1i.*(k-kc).*(L-l)./2)-1.)./(1i.*(k-kc)) +Ar(9).*Ccr.*(exp(1i.*(k+kc).*(L-l)./2)-1.)./(1i.*(k+kc)) +Ycc.*Ar(9).^2.*(exp(1i.*kc.*(L-l))-1.)./(2.*1i.*kc) + Ar(10).*Ccr.*(exp(1i.*(k-kc).*(L-l)./2)-1.)./(1i.*(k-kc)) -Ycc.*Ar(10).^2.*(exp(-1i.*kc.*(L-l))-1.)./(-2.*1i.*kc);
average_pvr=average_pvr./L^2;

%Compute average of the velocity: <v> for compressibility
% average_vc= 1./L^2*((Yct*Ac(1)/1i/kt)*(1 - exp(-1i*kt*L/2)) ...
% + (Yct*Ac(2,1)/1i/kt)*(1 -  exp(1i*kt*L/2)) ...
% + (Ctc/1i/k)*(1 - exp(-1i*k*L/2)) ...
% + (Yct*Ac(3,1)/1i/kt)*(exp(1i*kt*L/2) - 1) ...
% + (Yct*Ac(4,1)/1i/kt)*(exp(-1i*kt*L/2) - 1) ...
% + (Ctc/1i/k)*(exp(1i*k*L/2) - 1) ...
% + (Ycc*Ac(7,1)/1i/kc)*(1 - exp((-1/2)*1i*kc*(L-l))) ...
% + (Ycc*Ac(8,1)/1i/kc)*(1 - exp((1/2)*1i*kc*(L-l))) ...
% + (Ccc/1i/k)*(1 - exp((-1./2)*1i*k*(L-l))) ...
% + (Ycc*Ac(9,1)/1i/kc)*(exp((1/2)*1i*kc*(L-l)) - 1) ...
% + (Ycc.*Ac(10,1)./1i./kc).*(exp((-1./2).*1i.*kc.*(L-l)) - 1) ...
% + (Ccc./1i./k).*(exp((1./2).*1i.*k.*(L-l)) -1));

average_vc= Ctc.* ( exp(1i.*k.*L./2)-exp(-1i.*k.*L./2) ) ./1i./k +Ccc.*( exp(1i.*k.*(L-l)./2)-exp(-1i.*k.*(L-l)./2) ) ./1i./k +Yct.*Ac(1).*(1.-exp(-1i.*kt.*L./2))./1i./kt -Yct.*Ac(2).*(1.-exp(1i.*kt.*L./2))./(-1i)./kt +Ycc.*Ac(7).*(1.-exp(-1i.*kc.*(L-l)./2))./1i./kc -Ycc.*Ac(8).*(1.-exp(1i.*kc.*(L-l)./2))./(-1i)./kc + Yct.*Ac(3).*(exp(1i.*kt.*L./2)-1.)./1i./kt -Yct.*Ac(4).*(exp(-1i.*kt.*L./2)-1.)./(-1i)./kt +Ycc.*Ac(9).*(exp(1i.*kc.*(L-l)./2)-1.)./1i./kc -Ycc.*Ac(10).*(exp(-1i.*kc.*(L-l)./2)-1.)./(-1i)./kc;  
average_vc=average_vc./L^2;




%Compute <pv> for compressibility
% average_pvc= 1./(L^2)*(((Yct*Ac(1,1)^2)/(2*1i*kt))*(1 - exp(-1i*kt*L)) ...
% + ((Ac(1,1).*Ctc)./(1i.*(kt+k))).*(1 - exp((-1./2).*1i.*(kt+k).*L)) ...
% + ((Yct.*Ac(2,1).^2)./(2.*1i.*kt)).*(1 - exp(1i.*kt.*L)) ...
% + ((Ac(2,1).*Ctc)./(1i.*(k-kt))).*(1 - exp((-1./2).*1i.*(k-kt).*L)) ...
% + ((Yct.*Btc.*Ac(1,1))./(1i.*(kt+k))).*(1 - exp((-1./2).*1i.*(kt+k).*L)) ...
% - ((Yct.*Btc.*Ac(2,1))./(1i.*(k-kt))).*(1 - exp((-1./2).*1i.*(k-kt).*L)) ...
% + ((Btc.*Ctc)./(2.*1i.*k)).*(1 - exp(-1i.*k.*L)) ...
% + ((Yct.*Ac(3,1).^2)./(2.*1i.*kt)).*(exp(1i.*kt.*L) - 1) ...
% + ((Ac(3,1).*Ctc)./(1i.*(kt+k))).*(exp((1./2).*1i.*(kt+k).*L) -1) ...
% + ((Yct.*Ac(4,1).^2)./(2.*1i.*kt)).*(exp(-1i.*kt.*L) - 1) ...
% + ((Ac(4,1).*Ctc)./(1i.*(k-kt))).*(exp((1./2).*1i.*(k-kt).*L) - 1) ...
% + ((Yct.*Btc.*Ac(3,1))./(1i.*(kt+k))).*(exp((1./2).*1i.*(kt+k).*L) - 1) ...
% - ((Yct.*Btc.*Ac(4,1))./(1i.*(k-kt))).*(exp((1./2).*1i.*(k-kt).*L) - 1) ...
% + ((Btc.*Ctc)./(2.*1i.*k)).*(exp(1i.*k.*L) - 1) ...
% + ((Ycc.*Ac(7,1).^2)./(2.*1i.*kc)).*(1 - exp(-1i.*kc.*(L-l))) ...
% + ((Ac(7,1).*Ccc)./(1i.*(kc+k))).*(1 - exp((-1./2).*1i.*(kc+k).*(L-l))) ...
% + ((Ycc.*Ac(8,1).^2)./(2.*1i.*kc)).*(1 - exp(1i.*kc.*(L-l))) ...
% + ((Ac(8,1).*Ccc)./(1i.*(k-kc))).*(1 - exp((-1./2).*1i.*(k-kc).*(L-l))) ...
% + ((Ycc.*Bcc.*Ac(7,1))./(1i.*(kc+k))).*(1 - exp((1./2).*1i.*(kc+k).*(L-l))) ...
% - ((Ycc.*Bcc.*Ac(8,1))./(1i.*(k-kc))).*(1 - exp((-1./2).*1i.*(k-kc).*(L-l))) ...
% + ((Ycc.*Ac(9,1).^2)./(2.*1i.*kc)).*(exp(1i.*kc.*(L-l)) - 1) ...
% + ((Ac(9,1).*Ccc)./(1i.*(kc+k))).*(exp((1./2).*1i.*(kc+k).*(L-l)) - 1) ...
% + ((Ycc.*Ac(10,1).^2)./(2.*1i.*kc)).*(exp(-1i.*kc.*(L-l)) - 1) ...
% + ((Ac(10,1).*Ccc)./(1i.*(k-kc))).*(exp((1./2).*1i.*(k-kc).*(L-l)) - 1) ...
% + ((Ycc.*Bcc.*Ac(9,1))./(1i.*(kc+k))).*(exp((1./2).*1i.*(kc+k).*(L-l)) - 1) ...
% - ((Ycc.*Bcc.*Ac(10,1))./(1i.*(k-kc))).*(exp((1./2).*1i.*(k-kc).*(L-l)) - 1) ...
% + ((Bcc.*Ccc)./(2.*1i.*k)).*(exp(1i.*k.*(L-l)) - 1) ...
% + ((Bcc.*Ccc)./(2.*1i.*k)).*(1 - exp(-1i.*k.*(L-l))));


average_pvc=Btc.*Ctc.*(1.-exp(-1i.*k.*L))./(2.*1i.*k) +Ac(1).*Ctc.*(1.-exp(-1i.*(k+kt).*L./2))./(1i.*(k+kt)) +Ac(2).*Ctc.*(1.-exp(-1i.*(k-kt).*L./2))./(1i.*(k-kt)) +Btc.*Yct.*Ac(1).*(1.-exp(-1i.*(k+kt).*L./2))./(1i.*(k+kt)) +Yct.*Ac(1).^2.*(1.-exp(-1i.*kt.*L))./(2.*1i.*kt) -Btc.*Yct.*Ac(2).*(1.-exp(-1i*(k-kt).*L./2))./(1i.*(k-kt)) -Yct.*Ac(2).^2.*(1.-exp(1i.*kt.*L))./(-2.*1i.*kt);
average_pvc=average_pvc +Bcc.*Ccc.*(1.-exp(-1i.*k.*(L-l)))./(2.*1i.*k) +Ac(7).*Ccc.*(1.-exp(-1i.*(k+kc).*(L-l)./2))./(1i.*(k+kc)) +Ac(8).*Ccc.*(1.-exp(-1i.*(k-kc).*(L-l)./2))./(1i.*(k-kc)) +Bcc.*Ycc.*Ac(7).*(1.-exp(-1i.*(k+kc).*(L-l)./2))./(1i.*(k+kc)) +Ycc.*Ac(7).^2.*(1.-exp(-1i.*kc.*(L-l)))./(2.*1i.*kc) -Bcc.*Ycc.*Ac(8).*(1.-exp(-1i*(k-kc).*(L-l)./2))./(1i.*(k-kc)) -Ycc.*Ac(8).^2.*(1.-exp(1i.*kc.*(L-l)))./(-2.*1i.*kc);
average_pvc=average_pvc +Btc.*Ctc.*(exp(1i.*k.*L)-1.)./(2.*1i.*k) +Btc.*Yct.*Ac(3).*(exp(1i.*(k+kt).*L./2)-1.)./(1i.*(k+kt)) -Btc.*Yct.*Ac(4).*(exp(1i.*(k-kt).*L./2)-1.)./(1i.*(k-kt)) +Ac(3).*Ctc.*(exp(1i.*(k+kt).*L./2)-1.)./(1i.*(k+kt)) +Yct.*Ac(3).^2.*(exp(1i.*kt.*L)-1.)./(2.*1i.*kt) + Ac(4).*Ctc.*(exp(1i.*(k-kt).*L./2)-1.)./(1i.*(k-kt)) -Yct.*Ac(4).^2.*(exp(-1i.*kt.*L)-1.)./(-2.*1i.*kt);
average_pvc= average_pvc+Bcc.*Ccc.*(exp(1i.*k.*(L-l))-1.)./(2.*1i.*k) +Bcc.*Ycc.*Ac(9).*(exp(1i.*(k+kc).*(L-l)./2)-1.)./(1i.*(k+kc)) -Bcc.*Ycc.*Ac(10).*(exp(1i.*(k-kc).*(L-l)./2)-1.)./(1i.*(k-kc)) +Ac(9).*Ccc.*(exp(1i.*(k+kc).*(L-l)./2)-1.)./(1i.*(k+kc)) +Ycc.*Ac(9).^2.*(exp(1i.*kc.*(L-l))-1.)./(2.*1i.*kc) + Ac(10).*Ccc.*(exp(1i.*(k-kc).*(L-l)./2)-1.)./(1i.*(k-kc)) -Ycc.*Ac(10).^2.*(exp(-1i.*kc.*(L-l))-1.)./(-2.*1i.*kc);
average_pvc=average_pvc./L^2;




% average_bc= 1./(1i.*omega.*L.^2).*(2.*1i.*Ctc.*sin(kt.*L./2) ...
% - Yct.*Ac(1).*exp(-1i.*kt.*L./2)...
% + Yct.*Ac(2).*exp(1i.*kt.*L./2) ...
% + Yct.*Ac(3).*exp(1i.*kt.*L./2) ...
% - Yct.*Ac(4).*exp(-1i.*kt.*L./2));

average_bc=2.*1i.*Ctc.*sin(kt.*L./2) -Yct.*Ac(1,1).*exp(-1i.*kt.*L./2) +Yct.*Ac(2,1).*exp(1i.*kt.*L./2) +Yct.*Ac(3,1).*exp(1i.*kt.*L./2) -Yct.*Ac(4,1).*exp(-1i.*kt.*L./2);
average_bc=average_bc./(1i.*omega.*L^2);


% Compute H for density and compressibility
Hr= average_pvr./average_vr;
Hc= average_pvc./average_vc;


% Compute effective nonlocal density and compressibility

%rho= phi.*(f0-1i.*k.*Hr)./(-1i.*omega.*average_vr);
rho=phi.*(1.-1i.*k.*Hr)./(-1i.*omega.*average_vr);

chiinv=chi0inv.*(1.+Hc)./(average_bc.*(chi0inv./phi)+gamm);
chi=1./chiinv;

%rhoaf=porosite.*(1.-1i.*k.*Hrho)./(-1i.*omega.*vmoyrho);
































