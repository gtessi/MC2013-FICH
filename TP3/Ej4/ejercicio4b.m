clear all
close all
clc

syms x;

% condiciones de borde (Dirichlet)
izq=1;

% condiciones de borde (Neumann)
qder=0;

% cantidad de elementos
M=3;

% delta
dx=1/M;

% longitud
L=1;

% nodos
xnodos=0:dx:L;

% cantidad de nodos
n=length(xnodos);
%n=M+1;

Ke=zeros(2,2,M);
fe=zeros(2,M);

% calcula las matrices K de cada elemento
for e=1:M
    xi=xnodos(e);
    xj=xnodos(e+1);
    
    he=xj-xi;
    
    X=x-xi;
    
    N=[(he-X)/he X/he];
    
    dN=[-1/he 1/he];
    
    for ie=1:2
        for je=1:2
            integrando1=dN(ie)*dN(je);
            Ke(ie,je,e)=double(int(integrando1,x,xi,xj));
        end
    end
end

% matriz K global
K=zeros(n);

% vector f global
f=zeros(n,1);

% esambla la matriz global
for e=1:M
    % indices globales
    idx=[e e+1];
    
    % agrega la contribucion del elemento a la matriz global
    for ie=1:2
        ig=idx(ie);
        for je=1:2
            jg=idx(je);
            K(ig,jg)=K(ig,jg)+Ke(ie,je,e);
        end
    end
end

% la fuente interna Q es 1 si x<=L/2 y 0 si x>L/2
for ig=2:n-1
    xi=xnodos(ig);
    
    % si llega al punto medio de la malla, la integral cambia y reparte el
    % paso he a los nodos involucrados (cuando los nodos son pares)
    if abs(xi-L/2)<1e-3
        f(ig)=he/2;
        break;
    end
    
    % en los demas puntos menores al punto medio el paso es he
    f(ig)=he;
end

% condiciones de borde
% Dirichlet
K(1,:)=0;
K(1,1)=1;
f(1)=izq;

% Neumann
f(n)=-qder;

% resuelve el sistema
T=K\f;


