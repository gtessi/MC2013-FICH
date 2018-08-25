function N = N_tri(idx,x,y,nodosx,nodosy)
    % numeracion de los nodos a partir de idx
    i=idx;
    j=mod(idx,3)+1;
    k=mod(idx+1,3)+1;
    
    % valores de los nodos
    xi=nodosx(i);
    xj=nodosx(j);
    xk=nodosx(k);
    
    yi=nodosy(i);
    yj=nodosy(j);
    yk=nodosy(k);
    
    % coeficientes de la funcion de forma
    a=yk*xj-yj*xk;
    b=yj-yk;
    c=xk-xj;
    
    % area del triangulo
    A=det([1 xi yi; 1 xj yj; 1 xk yk]);
    
    % funcion de forma
    N=(a+b*x+c*y)/A;
end