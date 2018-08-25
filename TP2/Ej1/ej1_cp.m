function [K,f,a,error] = ej1_cp(M,L)
    % ------------------------------------------------
    % aproximacion de 1 + sin(pi*x/2) usando residuos ponderados, con
    % colocacion puntual
    % ------------------------------------------------
    % genera los puntos donde se evaluara el error
    x=linspace(0,L,(M+2)); % agrega los extremos a M
    x=x(2:(M+1)); % quita los extremos
    
    % define la matriz K y el vector f
    K=zeros(M);
    f=zeros(M,1);
    
    % calcula los valores de K y f
    for l=1:M
        % K
        for m=1:M
            K(l,m)=N(m,x(l),L);
        end
        % f
        f(l)=phi(x(l))-psi(x(l));
    end
    
    % resuelve el sistema
    a=K\f;
    
    % calcula el error
    x_e=linspace(0,L,1000);
    
    y_ex=phi(x_e);
    y_ap=psi(x_e);
    
    for m=1:M
        y_ap=y_ap+a(m)*N(m,x_e,L);
    end
    
    error=norm((y_ex-y_ap),2)/norm(y_ex,2);
    
    % grafica
    figure;
    hold on;
    plot(x_e,y_ex,'b');
    plot(x_e,y_ap,'r');
    title('Ejercicio 1 - Colocacion puntual');
    xlabel('x');
    ylabel('phi');
    legend('Exacta','Aproximada','Location','SouthEast');
end