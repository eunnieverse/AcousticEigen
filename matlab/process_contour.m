%---------------------------------------------------------------------
%- Process 'NewtCyl_phi90_Th100_rootfinding_s.txt'
%- Yoonkyung Eunnie Lee 2015. 08. 25 
%---------------------------------------------------------------------

data = dlmread('NewtCyl_phi90_Th100_rootfinding_s.txt'); 

klpi = data(1,1); 
disp (sprintf('at klpi = %g', klpi)); 
k1 = 36711.468838+(6592.761443i);
k2 = 525899.608420+(821968.653170i); 
disp (sprintf('answer k1 = %s, k2 = %s',num2str(k1,'%g'), num2str(k2,'%g'))); 
for ii= 1: size(data,1) % for each row of data 
    n = data(ii,2); 
    s0 = data(ii,3) + 1i * data(ii,4); 
    s1 = data(ii,5) + 1i * data(ii,6); 
    s2 = data(ii,7) + 1i * data(ii,8); 
        
    k2 = (-s1 - sqrt(2 * s2 - s1^2))/2;
    k1 = s1 - k2; 
    disp (sprintf('contour pts: %4d, k1 = %s, k2 = %s',n,num2str(k1,'%g'), num2str(k2,'%g'))); 
end