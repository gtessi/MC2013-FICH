function phi = diffin1DNeumann2(dx,L,ql,Tr,k,c,Q)
    % ------------------------------------------------
    % Metodo de Diferencias Finitas en 1D con condiciones Dirichlet
    % ------------------------------------------------
    % dx = diferencial x
    % L = longitud de la barra
    % ql = flujo del borde izquierdo
    % Tr = condicion del borde derecho
    % k = coeficiente de conductividad termica
    % c = coeficiente de reaccion
    % Q = fuente de calor
    % ------------------------------------------------
    %    d^2T
    % k*------ + Q - c*T = 0
    %    dx^2
    % ------------------------------------------------
    
    % nodos de la malla
    x=0:dx:L;
    
    % cantidad de nodos de la malla
    n=length(x)-2;
    
    % define los nodos centrales
    nodo=[1 -(2+(c/k)*dx.^2) 1];
    
    % define la matriz K
    K=zeros(n+1,n+1);
    f=ones(n+1,1)*(-Q/k)*dx.^2;
    
    % arma la matriz y el vector de terminos independientes
    % primer nodo
    K(1,(1:2))=[-(1+(c/(2*k))*dx.^2) 1];
    f(1)=(f(1)-(ql/k)*2*dx)/2;
    
    % nodos intermedios
    for l=2:n
        K(l,(l-1:l+1))=nodo;
    end
    
    % ultimo nodo
    K(n+1,(n:n+1))=nodo(1:2);
    f(n+1)=f(n+1)-Tr;
    
    % resuelve el sistema
    phi=K\f;
    
    % agrega los nodos de la frontera
    phi=[phi; Tr];
    
    % grafica la solucion aproximada
    figure;
    plot(x,phi);
    xlabel('x');
    ylabel('T');
end