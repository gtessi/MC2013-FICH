clear all
close all
clc

syms Xi;

% condiciones de borde (Dirichlet)
izq=0;
der=0;

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
    
    N=[(1-Xi)/2 (1+Xi)/2];
    
    dN=[-1/he 1/he];
    
    for ie=1:2
        for je=1:2
            integrando1=(dN(ie)*(2/he)*dN(je)*(2/he)+N(ie)*N(je))*(he/2);
            Ke(ie,je,e)=double(int(integrando1,Xi,-1,1));
        end
        %integrando2=N(ie)*(he*Xi/2)*(he/2);
        integrando2=N(ie)*(xi*N(1)+xj*N(2))*(he/2);
        fe(ie,e)=double(int(integrando2,Xi,-1,1));
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
        f(ig)=f(ig)+fe(ie,e);
    end
end

% condiciones de borde
% Dirichlet
K(1,:)=0;
K(1,1)=1;
f(1)=izq;

K(n,:)=0;
K(n,n)=1;
f(n)=der;

% resuelve el sistema
phi=K\f;



