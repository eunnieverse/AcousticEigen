% Find the macroscopic wave number - Helmholtz resonator
% Using Zwikker-Konten model for slits
clc;
clear;



% Air properties
cstphys3;
%cstphys3D;
%cstphys3_water;
% c0=343; 
% rho0=1.2;
% eta=1.84e-05;
% %eta=1.84e-20;
% p0=1.0132e+5;
% kappa=0.0262;
% T0=293;
% gamm=1.4;
% cp=1005;
% nup=kappa./(rho0.*cp);
% %Pr=eta./(rho0.*nup);
% Pr=0.71;


%wn=pi*wn^2/4;


Vc=wc*dc*hc;
%ht=ht*dt;
omegah=c0.*sqrt((pi*wn^2/4)./(hn.*Vc))              % resonance frequency 
omegap=sqrt( (Vc+ht*dt.*Lx)./(ht*dt.*Lx) ).*omegah

klpip=omegap.*Lx./c0./pi;

klpih=omegah.*Lx./(c0.*pi);


klpimin=0.00001;
klpimax=0.5;
pasklpi=(klpimax-klpimin)/10000;
pasomega=pasklpi.*c0.*pi./Lx;
omegamin=klpimin.*c0.*pi./Lx;
omegamax=klpimax.*c0.*pi./Lx;

%klpi=[klpimin: pasklpi: klpimax];

omega1=[omegamin: pasomega: omegap];
omega2=[omegap+pasomega: pasomega: omegamax];

klpi1=omega1.*Lx./c0./pi;
klpi2=omega2.*Lx./c0./pi;
klpi=[klpi1,klpi2];

omega=klpi.*c0.*pi./Lx;                  % angular frequency (Hz) 




c1w= c0./sqrt( 1-(Vc.*omegah.^2)./(ht*dt.*Lx.*(omega1.^2-omegah.^2)) );
c2w= c0./sqrt( 1-(Vc.*omegah.^2)./(ht*dt.*Lx.*(omega2.^2-omegah.^2)) );
cw=[c1w,c2w];         % phase velocity in ideal fluid computed by the "formula"

PLOT_k=plot(real(omega./cw), omega/(2*pi)*1e-3, 'b-'); hold on;
plot(imag(omega./cw), omega/(2*pi)*1e-3, 'r-');

 
 
 
 
