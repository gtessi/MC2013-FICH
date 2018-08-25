function K = Ensambler(icone,Ke)
    % ------------------------------------------------------
    % ENTRADA
    % icone = matriz que relaciona elementos con nodos
    % Ke = vector de matrices elementales
    % ------------------------------------------------------
    % SALIDA
    % K = matriz global ensamblada
    % ------------------------------------------------------
    
    % cantidad de elementos
    ne=size(icone,1);
    
    % cantidad de nodos
    nn=max(max(icone))*2;
    
    % K global
    K=zeros(nn);
    
    % recorre los elementos
    for e=1:ne
        % nodos del elemento
        nodos=icone(e,:);
        
        % matriz del elemento
        elemento=Ke(:,:,e);
        
        K(nodos,nodos)=K(nodos,nodos)+elemento;
    end
end