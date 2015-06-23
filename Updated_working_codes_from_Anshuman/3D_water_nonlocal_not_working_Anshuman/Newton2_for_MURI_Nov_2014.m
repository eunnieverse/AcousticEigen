% Compute Nonlocal wavenumbers - Array of Helmholtz resonators
% Using Zwikker-Konten model for slits
% Taking into account Zwikker-Konten waves in the cavity-- Vertically
% isOpen = matlabpool('size') > 0;
% if ~isOpen 
%     matlabpool open
% end



% Geometry
clear all
clc;
%%%%%clearvars -global


% Physical constants
%constantes
%cstphys3_water;
cstphys3;
err_bound=1e-5;
%nu
%rho0
%chi0inv
% Compute nonlocal wavenumber for each frequency by Newton-Raphson

klpimin= 0.00001;
%klpimin= 0.001;
klpimax= 0.5;
nbptklpi=100;
klpi=linspace(klpimin,klpimax,nbptklpi);



omega=klpi*pi*c0 /Lx;
chiNL=0*omega;
rhoNL=0*omega;
KNL=0*omega;
kNL=0*omega;


% Block wavenumber
[q]=Wavenumber2(omega,nbptklpi);


% dk
%eps= (1.e-11 +0.*1i).*q;
eps=1.e-6 +0.*1i;

qin=zeros(nbptklpi,1); 
qout=zeros(nbptklpi,1); 

%commented this part for testing the Bloch code
tag=0;
Ftoplot=zeros(length(klpi));
for nn=1:nbptklpi
%   
nn
%     if nn==1
%        qin(nn)=q(nn).*(1.+0.01);
%     else
%        qin(nn)=qin(nn-1);
%      end
   qin(nn)= q(nn).*(1.+0*0.05);   % Initial value of wavenumber for Newton-Raphson 
    err=1+1i;
    %count=1;
    while abs(err)>err_bound
    %[valrhoAf, valchiAf]=qattf(omega(nn),qin(nn));
    
    %%%%%%%%%%These two lines were removed on 1st Aug, 2014
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho, chi]=Density2(omega(nn),qin(nn));
[rho1, chi1]=Density2(omega(nn),qin(nn).*(1.+eps));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %     Dchi=valchiAf-chi;
%     Drho=valrhoAf-rho;
    F= rho.*chi.*omega(nn).^2 - qin(nn).^2;                                % Nonlocal dispersion equation: F=0
    F1= rho1.*chi1.*omega(nn).^2 - (qin(nn).*(1.+eps)).^2;                 % F at k+dk
    Fp=(F1-F)./(qin(nn).*eps);                                                     % Derivative of F     
    
    qout(nn)=qin(nn)-F./Fp;
    %Fplot(count)=F;
    %Kplot(count)=qin(nn);
    %err=(qout(nn)-qin(nn))./qin(nn);
    err=abs(F);
    qin(nn)=qout(nn);
    %count=count+1;
    %if count>1000
     %   qin(nn)=q(nn).*(1.+0.1);
      %  tag=tag+1;
       % count=1;
        %Fplot=[];
        %Kplot=[];
    %end
    end
%     if imag(qout(nn))<0
%         qout(nn)=-qout(nn);
%     end
%     Ftoplot(nn)=F;
end

 kNL= qout;
lambdaNL= 2*pi./real(kNL)/Lx; 
lambda0=2*pi./omega.*(c0/Lx); 







%hold on


for nn=1:nbptklpi
    valomega=omega(nn);
    valknl=kNL(nn);
    
[rho, chi]=Density2(valomega,valknl);

rhoNL(nn)=rho;
chiNL(nn)=chi;
KNL(nn)= 1./chi;
end


figure;
PLOT_k=plot(klpi,real(q), 'b-');hold on;
plot(klpi, imag(q),'r-')
xlim([0 0.5])
PLOT_kNL=plot(klpi,real(kNL), 'bo');hold on;
plot(klpi, imag(kNL),'ro')

hleg = legend( [PLOT_k, PLOT_kNL], 'Bloch-wave model', 'Nonlocal model');
xlabel(' k_0 L/\pi')
 ylabel('k (m^{-1})')
  hold off;

 % PLOT rho
figure;
% Plot the model results
 PLOT_rho_nonlocal = plot(klpi(10:end),real(rhoNL(10:end)),'b-');
 hold on
 plot(klpi(10:end),imag(rhoNL(10:end)),'r-')

xlabel('k_0 L/\pi' )
ylabel('\rho (kg m^{-3})' )


hold off;
figure;
% Plot the model results
PLOT_chi_nonlocal = plot(klpi(3:end),real(1./chiNL(3:end)),'b-');
hold on
plot(klpi(3:end),imag(1./chiNL(3:end)),'r-')

xlabel('k_0 L/\pi' )
ylabel('K (Pa)' )

