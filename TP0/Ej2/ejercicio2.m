close all
clear all
clc
% ------------------------------------------------
% datos
% longitud de la malla
L=5;
% paso de la malla
dx=0.5;
% nodos de la malla
x=0:dx:L;

% factores
A=(-1-exp(-L))/(exp(L)+exp(-L));
B=(1-exp(L))/(exp(L)+exp(-L));

% solucion analitica
y=exp(x)*A+exp(-x)*B+1;

% grafica la solucion analitica
figure;
plot(x,y);
title('Ejercicio 2 - Soluciones');
xlabel('x');
ylabel('T');
% ------------------------------------------------
% valores de nodos en la frontera
Tl=0;
qr=1; % condicion neumann

% cantidad de nodos de la malla
n=length(x)-2;

% define los nodos centrales
nodo=[1 -(2+dx.^2) 1];

% define la matriz K
K=zeros(n+1,n+1);
f=zeros(n+1,1);

% arma la matriz y el vector de terminos independientes
% primer nodo
K(1,(1:2))=nodo(2:3);
f(1)=-(dx.^2+Tl);

% nodos intermedios
for k=2:n
    K(k,(k-1:k+1))=nodo;
    f(k)=-dx.^2;
end

% ultimo nodo
K(n+1,(n:n+1))=[1 -(2+dx.^2)/2];
f(n+1)=(-dx.^2+2*dx*qr)/2;

% resuelve el sistema
phi=K\f;

% agrega los nodos de la frontera
phi=[Tl; phi];

% grafica la solucion aproximada
hold on;
plot(x,phi,'r');
legend('Solucion analitica','Solucion aproximada','location','SouthEast');
% ------------------------------------------------