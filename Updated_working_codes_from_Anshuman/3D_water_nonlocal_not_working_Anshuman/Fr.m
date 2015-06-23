function [ out ] = Fr(mu,omega,a,b)
%# for rectangular tube: main tube and cavity; dimensions: 2a and 2b
err=1e-4;    

out = 0*omega;
for l= 1 : length(omega)
    om=omega(l);
    
coef=-4j*om/mu/a^2/b^2;

    m=0;
    sum2=0;
    errm=1;
    while errm>err
        n=0 ;          
        sum1=0;
        errn=1;
        while errn>err
            add1=1/alpha(m,a)^2/beta(n,b)^2/(alpha(m,a)^2+beta(n,b)^2-1j*om/mu);
            %errn=abs(add1);
            sum1=sum1+add1;
            errn=abs(add1./sum1);
            n=n+1;
        end
        sum2=sum2+sum1;
        errm=abs(sum1./sum2);
        m=m+1;
    end
    out(l)= sum2*coef;

end

