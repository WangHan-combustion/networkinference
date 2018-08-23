function [p]=lognormal(x,mu,var)

pi=22.0/7.0;

p=-0.5*log(2*pi*var)-((x-mu).^2)/(2*var);