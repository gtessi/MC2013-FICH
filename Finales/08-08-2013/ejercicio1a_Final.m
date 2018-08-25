close all
clear all
clc
% ------------------------------------------------
% datos
% paso
h=0.25;

% nodos de la malla no uniforme
x=0:h:1;

% cantidad de nodos
n=length(x);

% condiciones de borde
izq=1;
der=2;

% matriz global
K=zeros(n,n);

% vector de constantes
f=zeros(n,1);

% calcula los nodos en la frontera
% izquierdo
K(1,1)=1;
f(1)=izq;

% derecho
K(n,n)=1;
f(n)=der;

% calcula los nodos internos de la malla
for k=2:n-1
    hp=x(k+1)-x(k);
    hn=x(k)-x(k-1);
    K(k,(k-1:k+1))=[2/(hn*(hp+hn)) (1-2/(hp*hn)) 2/(hp*(hp+hn))];
    f(k)=10*x(k)*(x(k)+1);
end

% resuelve el sistema
phi=K\f;

% grafica la solucion aproximada
figure;
plot(x,phi)
xlabel('x');
ylabel('T');
% ------------------------------------------------