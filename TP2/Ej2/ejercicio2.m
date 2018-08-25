close all
clear all
clc
% ------------------------------------------------
% datos
% largo de la barra
L=1;

% cantidad de nodos
M=1;

% variable simbolica
syms x;

% define la matriz K y el vector f
K=zeros(M);
f=zeros(M,1);

% calcula los valores de K y f
for l=1:M
    % K
    for m=1:M
        K(l,m)=double(int((DN(l,x)*DN(m,x)-N(l,x)*N(m,x)),x,0,L)+N(l,L)*N(m,L));
    end
    % f
    f(l)=double(int(N(l,x),0,L));
end

% resuelve el sistema
a=K\f;

% calcula el error
x_e=linspace(0,L,1000);

y_ex=phi(x_e);
y_ap=psi(x_e);

for m=1:M
    y_ap=y_ap+a(m)*N(m,x_e);
end

error=norm((y_ex-y_ap),2)/norm(y_ex,2);

% grafica
figure;
hold on;
plot(x_e,y_ex,'b');
plot(x_e,y_ap,'r');
title('Ejercicio 2 - Galerkin');
xlabel('x');
ylabel('phi');
legend('Exacta','Aproximada','Location','SouthEast');
% ------------------------------------------------