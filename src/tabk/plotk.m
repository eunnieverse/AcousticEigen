function cfig=plotk(datfilename)
%---------------------------------------------------------------------
%- plots klpi1, klpi2 data for given index list 
%- Yoonkyung Eunnie Lee 2015. 08. 30
%---------------------------------------------------------------------
% choose plot color 
clist = 'bkrgcm';
sty = sprintf('%s-o',clist(1)); 'k-o';
%---------------------------------------------------------------------
% Load Data 
data = dlmread(datfilename); 
index  = data(:,1); 
klpi = data(:,2); 
k =data(:,3)+ 1i.*data(:,4);
k = k/100000;
%sampleii = 17:70; 
sampleii=1:length(data);
%---------------------------------------------------------------------
cfig = figure(); 
set(cfig,'Position',[20 20 20+920 20+460]);
subplot(2,2,[1,3]);
 plot(real(k(sampleii)),imag(k(sampleii)),sty);
 xlim([-2 18]);  ylim([-2 18]); axis square; 
 xlabel('Re(k) [10^5 rad/m]'); ylabel('Im(k)  [10^5 rad/m]'); 
 title(datfilename,'Interpreter','none');
 
subplot(2,2,2); 
plot(klpi(sampleii),real(k(sampleii)),sty); hold on; 
axis tight;
ylabel('Re(k)'); 
subplot(2,2,4); 
plot(klpi(sampleii),imag(k(sampleii)),sty); hold on; 
 axis tight;
 ylabel('Im(k)'); 
 xlabel('kL/\pi, nondimentionalized frequency'); 
end