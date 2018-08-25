clear all
close all
clc

% condiciones de borde (Dirichlet)
izq=0;
der=1;

% cantidad de elementos
M=3;

% delta
dx=1/M;

% longitud
L=1;

% nodos
x=0:dx:L;

% cantidad de nodos
n=length(x);
%n=M+1;

% matriz K global
K=zeros(n);

% vector f global
f=zeros(n,1);

% esambla la matriz global
for k=1:M
    % calcula he
    he=x(k+1)-x(k);
    
    % elemento K
    Ke=[1/he+he/3 -1/he+he/6; -1/he+he/6 1/he+he/3];
    
    % dimension del elemento
    ne=size(Ke,1);
    
    % agrega la contribucion del elemento a la matriz global
    K(k:k+ne-1,k:k+ne-1)=K(k:k+ne-1,k:k+ne-1)+Ke;
end

% reduce el sistema sin la primera ni la ultima fila y columna
Kred=K(2:M,2:M);

fred=zeros(ne,1);
fred(1)=-K(2,1)*izq;
fred(M-1)=-K(M+1,M)*der;

% calcula las incognitas
phi=Kred\fred;

% agerga los bordes
phi=[0; phi; 1];

% agrega las condiciones de borde (Dirichlet)
f(1)=-(K(1,1:ne)*phi(1:ne));
f(M+1)=(K(M+1,M+1-ne:M+1)*phi(M+1-ne:M+1));


