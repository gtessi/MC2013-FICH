clear all

F3x=100000;
F3y=-50000;

E=20*10^9;

nu=0.3;

D=(E/(1-nu.^2)).*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

h=1;

A=(h.^2)/2;

xnod=[0 h; h 2*h; h h; 0 0; 2*h 2*h; 2*h h; h 0];

% la numeracion no se hace asi, para la coordenada x es numero de nodo * 2
% - 1 y para la coordenada y es numero de nodo * 2
icone=[1 2 5 6 3 4;
       7 8 5 6 1 2;
       5 6 9 10 3 4;
       7 8 13 14 5 6;
       5 6 11 12 9 10;
       13 14 11 12 5 6];

% parte elemental tipo 1
B11=[-1 0; 0 0; 0 -1];
B12=[1 0; 0 -1; -1 1];
B13=[0 0; 0 1; 1 0];
B1=[B11 B12 B13];

% parte elemental tipo 2
B21=[0 0; 0 -1; -1 0];
B22=[1 0; 0 0; 0 1];
B23=[-1 0; 0 1; 1 -1];
B2=[B21 B22 B23];

% matrices elementales
vec1=[1 4 5];
vec2=[2 3 6];

Ke=zeros(6,6,6);

for e=1:3
    % elementos tipo 1
    Ke(:,:,vec1(e))=B1'*D*B1*A;
    
    % elementos tipo 2
    Ke(:,:,vec2(e))=B2'*D*B2*A;
end

ne=size(icone,1);
    
% cantidad de nodos
nn=max(max(icone));

% K global
Kg=zeros(nn);

% recorre los elementos
for e=1:ne
    % nodos del elemento
    nodos=icone(e,:);
    
    % matriz del elemento
    elemento=Ke(:,:,e);
    
    % ensambla la matriz
    Kg(nodos,nodos)=Kg(nodos,nodos)+elemento;
end

% arma el vector de cargas
f=[0 0 0 0 F3x F3y 0 0 0 0 0 0 0 0]';

% resuelve el sistema
u=Kg\f;
