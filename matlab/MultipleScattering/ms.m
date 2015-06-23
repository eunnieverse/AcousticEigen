load('Navid_Phi90_DATA_PTS100.mat'); 
c1 = transpose(c1);
k1 = transpose(K1); 
clear K1 c2 K2 c3 K3 Vec_kalpi ;
data =dlmread('tabk_to_ms_c.txt');
klpi_ms=data(:,1);
c2= data(:,2) + i*data(:,3);
data =dlmread('tabk_to_ms_k.txt');
k2= data(:,2) + i*data(:,3);
clear data
save('ms.mat','c1','c2','k1','k2','klpi_ms');