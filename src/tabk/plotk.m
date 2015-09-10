function cfig = plotk()
%---------------------------------------------------------------------
%- plots klpi1, klpi2 data for given index list 
%- Yoonkyung Eunnie Lee 2015. 08. 30
%---------------------------------------------------------------------
tabk1 = dlmread('NewtCyl_phi90_tabk_Mode1.dat'); 
klpi  = tabk1(:,4); 
tabk1 = tabk1(:,2)+1i*tabk1(:,3);
tabk2 = dlmread('NewtCyl_phi90_tabk_Mode2.dat'); 
tabk2 = tabk2(:,2)+1i*tabk2(:,3);
sampleii = 17:70;
%---------------------------------------------------------------------
cfig=figure();
 plot(real(tabk1(sampleii)),imag(tabk1(sampleii)),'b-*'); hold on; 
 plot(real(tabk2(sampleii)),imag(tabk2(sampleii)),'r-*'); 
 legend('k1','k2'); axis equal; axis tight;
 xlabel('Re(k)'); ylabel('Im(k)'); 
%---------------------------------------------------------------------
end