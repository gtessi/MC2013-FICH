function F = Fsambler(icone,fe)
    % ------------------------------------------------------
    % ENTRADA
    % icone = matriz que relaciona elementos con nodos
    % fe = vector de vectores independientes
    % ------------------------------------------------------
    % SALIDA
    % F = vector global ensamblado
    % ------------------------------------------------------
    
    % cantidad de elementos
    ne=size(icone,1);
    
    % cantidad de nodos
    nn=max(max(icone));
    
    % F global
    F=zeros(nn,1);
    
    % recorre los elementos
    for e=1:ne
        % nodos del elemento
        nodos=icone(e,:);
        
        % matriz del elemento
        elemento=fe(:,:,e);
        
        F(nodos,1)=F(nodos,1)+elemento;
    end
end