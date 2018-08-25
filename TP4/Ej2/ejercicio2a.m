clear all
close all
clc

syms x y;

% ---------------------------------------
% constantes
Kelvin=273.15;

% datos
% cantidad de elementos
E=2;

% ancho y alto de la placa
Lx=5;
Ly=5;
% deltas
dx=5;
dy=5;
% nodos de la placa
nodosx=0:dx:Lx;
nodosy=0:dy:Ly;

% espesor
t=1;
% conductividad termica del material
k=2;
D=k*eye(2);

% condiciones de borde
% frontera 1 (inferior)
% *** aislada ***
% frontera 2 (lateral derecha)
T2_0=100+Kelvin;
% frontera 3 (superior)
h=1.2;
T3_amb=30+Kelvin;
% frontera 4 (lateral izquierda)
q4=2;

% fuentes
% elemento 1
Qp=5; % fuente puntual en P=(1,1)
% elemento 2
Q=1.2;
% ---------------------------------------

% elemento 1

he=5;

A=(he.^2)/2;

B=[-he he 0; -he 0 he]./(2*A);

N=[he.^2-he*x-he*y 0+he*x+0 0+0+he*y]./(2*A);

% matriz K
% integral sobre el dominio
Ke(:,:,1)=B'*D*B*t*A;

% vector f
fe(:,1)=-q4.*int(subs(N.',x,0),y,0,5)+Qp*subs(subs(N.',x,1),y,1);


% elemento 2

he=5;

A=(he.^2)/2;

B=[0 he -he; -he he 0]./(2*A);

N=[he.^2+0-he*y -he.^2+he*x+he*y he.^2-he*x+0]./(2*A);

% matriz K
% integral sobre el dominio
Ke(:,:,2)=B'*D*B*t*A+h*int(subs(N.',y,5)*subs(N,y,5),x,0,5);

% vector f
fe(:,2)=Q*int(int(N.',x,0,5),y,0,5)%-h*T3_amb*int(subs(N.',y,5),x,0,5);


% *** armarlo a pata y comparar con el script calor2D.m
