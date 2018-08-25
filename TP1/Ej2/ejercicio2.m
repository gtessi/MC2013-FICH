close all
clear all
clc
% ------------------------------------------------
% datos
% largo de la barra
L=1;

% paso de la malla
dx=0.1;

% paso al cuadrado
dx2=dx.^2;

% nodos de la malla
x=0:dx:L;

% cantidad de nodos de la malla
n=length(x);

% define la matriz K y el vector f
K=zeros(n,n);
f=zeros(n,1);

% estructura para cada nodo
nodo=[1 -2+dx2 1];

% crea la matriz y el vector
for k=2:n-1
    K(k,(k-1:k+1))=nodo;
    f(k)=4*dx2*x(k)*(x(k)-1);
end

% corrige los nodos de la frontera porque hay flujos en los extremos
K(1,1:2)=[(-2+dx2)/2 1];
f(1)=2*dx2*x(k)*(x(k)-1)-10*dx;

K(n,(n-1:n))=[1 (-2+dx2)/2];
f(n)=2*dx2*x(k)*(x(k)-1);

% resuelve el sistema
phi=K\f;

% grafica la solucion aproximada
figure;
plot(x,phi);
xlabel('x');
ylabel('T');
% ------------------------------------------------