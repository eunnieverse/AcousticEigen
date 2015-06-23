function [ out ] =  Fc(mu,om)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%def Fc(mu,om):  # for circular tube with the radius rw
cstphys3;
    out= 1 - 2./rw*(1j*om/mu).^(-0.5) .*besselj(1,rw*(1j*om./mu).^(0.5))./besselj(0,rw*(1j*om./mu).^(0.5));
end

