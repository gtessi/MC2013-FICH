close all
clear all
clc
% ------------------------------------------------
% datos
% limites de la placa
% x
inf_x=-1;
sup_x=1;
% y
inf_y=-1;
sup_y=1;

% deltas
%dx=0.25;
%dy=0.25;
dx=0.05;
dy=0.05;

% cantidad de nodos
M=2;

% variables simbolicas
syms x y;

% define la matriz K y el vector f
K=zeros(M);
f=zeros(M,1);

% calcula los valores de K y f
for l=1:M
    % K
    for m=1:M
        K(l,m)=double(int(int(N(l,x,y)*(D2N(m,x,y)+D2N(m,x,y)),x,inf_x,sup_x),y,inf_y,sup_y));
    end
    % f
    f(l)=4*double(int(int(N(l,x,y),x,inf_x,sup_x),y,inf_y,sup_y));
end

% resuelve el sistema
a=K\f;

% define los intervalos de valores para x e y
x=inf_x:dx:sup_x;
y=inf_y:dy:sup_y;

% genera la malla
[X,Y]=meshgrid(x,y);

% calcula la funcion aproximada
T=psi(X,Y);

for m=1:M
    T=T+a(m)*N(m,X,Y);
end

% grafica
figure;
surf(X,Y,T);
title('Ejercicio 3 - Galerkin');
xlabel('x');
ylabel('y');
zlabel('T');
% ------------------------------------------------