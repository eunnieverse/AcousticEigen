% Find the macroscopic wave number - Helmholtz resonator
% Using Zwikker-Konten model for slits
% Taking into account Zwikker-Konten waves in the cavity-- Vertically



function [q]=Wavenumber2(omega, nbptklpi)
%constantes
cstphys3

st=sqrt( omega.*rho0.*((ht./2).^2)./eta );                   % define parameter
sn=sqrt( omega.*rho0.*((wn./2).^2)./eta );                   % define parameter
sc=sqrt( omega.*rho0.*((hc./2).^2)./eta );                   % define parameter


rhot = rho0./( 1-tanh(st.*sqrt(-1i))./(st.*sqrt(-1i)) );         % effective density of the main tube (slit)          
rhon=  rho0./( 1-tanh(sn.*sqrt(-1i))./(sn.*sqrt(-1i)) );         % effective density of the neck (slit)
rhoc= rho0./( 1-tanh(sc.*sqrt(-1i))./(sc.*sqrt(-1i)) );          % effective density of the cavity (slit)

Kc=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*sc.*sqrt(-1i))./(sqrt(Pr).*sc.*sqrt(-1i)) );              % effective modulus of the cavity (slit)
Kt=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*st.*sqrt(-1i))./(sqrt(Pr).*st.*sqrt(-1i)) );               % effective modulus of the tube (slit)
Kn=  gamm.*p0./( 1 + (gamm-1).*tanh(sqrt(Pr).*sn.*sqrt(-1i))./(sqrt(Pr).*sn.*sqrt(-1i)) );              % effective modulus of the neck (slit)
  
%chit=1./Kt;
%chin=1./Kn;
%chic=1./Kc;

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

% for horizental waves in the cavity

Y6 = -2.*Ycc.* ( 1-exp(-1i.*kc.*(L-l)) )./ (1 + exp(-1i.*kc.*(L-l)) );
Yr= (-1i.*Ycn.*sin(kn.*l)+Y6.*cos(kn.*l))./(cos(kn.*l)-1i.*Y6.*sin(kn.*l)./Ycn);


% solve a second degree equation to get exp(iqL), with coefficients A, B and C



% Vertical waves
%A = Zr.*Zct;
%B = 1i.*(Zct.^2).*sin(kt.*L) - 2.*Zr.*Zct.*cos(kt.*L);
%C = Zr.*Zct;


% Horizontal waves
A = 1.;
B = 1i.*(Yr./Yct).*sin(kt.*L)-2.*cos(kt.*L);
C = 1;



q1= -(1i./L).*log((-B + sqrt(B.^2-4.*A.*C))./(2.*A));           % First solution
q2 = -(1i./L).*log((-B - sqrt(B.^2-4.*A.*C))./(2.*A));           % Second solution

% We look for propagating modes with wavenumber q, i.e., the modes which have positive imaginary part  

for nn=1: nbptklpi

if imag(q1(nn))>0
    q(nn)=q1(nn);  
else
    q(nn)=q2(nn);
end

end

    


%semilogx(kLpiNL,imag(cNL),'ob')
%legend('Multiple scattering','FEM local', 'FEM nonlocal')
