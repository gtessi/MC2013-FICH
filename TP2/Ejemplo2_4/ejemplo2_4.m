clear all
close all
clc

M=3;

K=zeros(M);
f=zeros(M,1);

Linf=0;
Lsup=1;

syms x;

for l=1:M
    f(l)=N(l,1);
    for m=1:M
        K(l,m)=double(int(-N(l,x)*(D2N(m,x)-N(m,x))+N(l,0)*N(m,0)+N(l,1)*N(m,1),x,Linf,Lsup));
    end
end

a=K\f;