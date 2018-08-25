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
A=exp(-L)/(exp(L)-exp(-L));
B=exp(L)/(exp(-L)-exp(L));

% solucion analitica
y=exp(x)*A+exp(-x)*B+1;

% grafica la solucion analitica
figure;
plot(x,y);
title('Ejercicio 1 - Soluciones');
xlabel('x');
ylabel('T');
% ------------------------------------------------
% valores de nodos en la frontera
Tl=0;
Tr=1;

% cantidad de nodos de la malla
n=length(x)-2;

% define los nodos centrales
nodo=[1 -(2+dx.^2) 1];

% define la matriz K
K=zeros(n,n);
f=zeros(n,1);

% arma la matriz y el vector de terminos independientes
% primer nodo
K(1,(1:2))=nodo(2:3);
f(1)=-(dx.^2+Tl);

% nodos intermedios
for k=2:n-1
    K(k,(k-1:k+1))=nodo;
    f(k)=-dx.^2;
end

% ultimo nodo
K(n,(n-1:n))=nodo(1:2);
f(n)=-(dx.^2+Tr);

% resuelve el sistema
phi=K\f;

% agrega los nodos de la frontera
phi=[Tl; phi; Tr];

% grafica la solucion aproximada
hold on;
plot(x,phi,'r');
legend('Solucion analitica','Solucion aproximada','location','SouthEast');
% ------------------------------------------------
% flujo de calor
q=-exp(x)*A+exp(-x)*B;

% grafica el flujo de calor
figure;
plot(x,q);
title('Ejercicio 1 - Flujo de calor');
xlabel('x');
ylabel('T');
% ------------------------------------------------