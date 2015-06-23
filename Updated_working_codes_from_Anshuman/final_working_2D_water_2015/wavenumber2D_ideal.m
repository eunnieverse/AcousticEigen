% Find the macroscopic wave number - Helmholtz resonator
% Using Zwikker-Konten model for slits
clear all;
clc


% Geometry
cstphys3_water; % change this for air, if required
% L=1e-2;
% ht=0.2*L; % height of main tube (m)
% hn=0.15*L; % height of neck (m)
% wn=0.015*L; % width of neck (m)
% Lx=L; % length of cell /x (m)
% Ly=L; % length of cell /y (m)
% wc=L-hn; % width of cavity (m)
% hc=L-2*hn-ht; % height of cavity (m)
Vc= hc.*wc; %volume of the cavity

% Air properties
%cstphys3
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

omegah=c0.*sqrt(wn./(hn.*Vc));              % resonance frequency 
omegap=sqrt( (Vc+ht.*Lx)./(ht.*Lx) ).*omegah;

klpip=omegap.*Lx./c0./pi;

klpih=omegah.*Lx./(c0.*pi);


klpimin=0.002;
klpimax=2;
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



st=sqrt( omega.*rho0.*((ht./2).^2)./eta );                   % define parameter
sn=sqrt( omega.*rho0.*((wn./2).^2)./eta );                   % define parameter
sc=sqrt( omega.*rho0.*((wc./2).^2)./eta );                   % define parameter


rhot = rho0./( 1-tanh(st.*sqrt(-1i))./(st.*sqrt(-1i)) );         % effective density of the main tube (slit)          
rhon=   rho0./( 1-tanh(sn.*sqrt(-1i))./(sn.*sqrt(-1i)) );        % effective density of the neck (slit)

Kc=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*sc.*sqrt(-1i))./(sqrt(Pr).*sc.*sqrt(-1i)) );            % effective modulus of the cavity (slit)
Kt= gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*st.*sqrt(-1i))./(sqrt(Pr).*st.*sqrt(-1i)) );               % effective modulus of the tube (slit)
Kn=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*sn.*sqrt(-1i))./(sqrt(Pr).*sn.*sqrt(-1i)) );             % effective modulus of the neck (slit)

ct= sqrt(Kt./rhot);        % phace velocity in the tube
cn= sqrt(Kn./rhon);        % phace velocity in the neck

kt= omega./ct;              % wave number of the tube (slit)
kn= omega./cn;              % wave number of the neck (slit)

Zct = rhot.*ct./ht;         % characteristic impedance of the main tube
Zcn= rhon.*cn./wn;          % characteristic impedance of the neck
Zr = ( -cos(kn.*hn)+(omega.*Vc./Kc).*Zcn.*sin(kn.*hn) )./( 1i.*sin(kn.*hn)./Zcn...
    + (1i.*omega.*Vc./Kc).*cos(kn.*hn) );

% solve a second degree equation to get exp(iqL), with coefficients A, B and C

A = Zr.*Zct;
B = 1i.*(Zct.^2).*sin(kt.*Lx) - 2.*Zr.*Zct.*cos(kt.*Lx);
C = Zr.*Zct;




q1 = -(1i./Lx).*log((-B + sqrt(B.^2-4.*A.*C))./(2.*A));           % First solution
q2 = -(1i./Lx).*log((-B - sqrt(B.^2-4.*A.*C))./(2.*A));           % Second solution

% We look for propagating modes with wavenumber q, i.e., the modes which have positive imaginary part  

for nn=1:max(size(omega));

if imag(q1(nn))>0
    q(nn)=q1(nn);  
else
    q(nn)=q2(nn);
end

end

    
c=omega./q;     % Phase velocity of the propagating modes


c1w= c0./sqrt( 1-(Vc.*omegah.^2)./(ht.*Lx.*(omega1.^2-omegah.^2)) );
c2w= c0./sqrt( 1-(Vc.*omegah.^2)./(ht.*Lx.*(omega2.^2-omegah.^2)) );
cw=[c1w,c2w];         % phase velocity in ideal fluid computed by the "formula"

qlpi=real(q).*Lx./pi;      % compute qL/pi to see if it is smaller than 0.5
figure;
xlim([0 klpimax])
PLOT_k=plot(klpi,real(omega./cw), 'b-');hold on;
xlim([0 klpimax])
plot(klpi, imag(omega./cw),'r-')


xlabel(' k_0 L/\pi' )
 ylabel('k (m^{-1})' )

hold off;


 
 
