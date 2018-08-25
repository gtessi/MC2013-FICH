function F = movie(metodo,Lx,Ly,Lt,dx,dy,dt)
    % ------------------------------------------------
    % funcion que crea una secuencia de la evolucion de las graficas
    % ------------------------------------------------
    if (lower(metodo)=='e')
        [phi,Nx,Ny,X,Y]=FE(Lx,Ly,Lt,dx,dy,dt);
    elseif (lower(metodo)=='i')
        [phi,Nx,Ny,X,Y]=BE(Lx,Ly,Lt,dx,dy,dt);
    elseif (lower(metodo)=='c')
        [phi,Nx,Ny,X,Y]=CN(Lx,Ly,Lt,dx,dy,dt);
    end
    
    nt=size(phi,2);
    
    F(nt)=struct('cdata',[],'colormap',[]);
    figure;
    
    for t=1:nt
        m_phi=reshape(phi(:,t),Nx+1,Ny+1);
        surf(X,Y,m_phi);

        axis([0,Lx,0,Ly,0,50]);
        grid on;
        set(gca,'xtick',0:dx:Lx);
        set(gca,'ytick',0:dy:Ly);
        set(gca,'ztick',0:10:50);
        xlabel('x');
        ylabel('y');
        zlabel('u');
        
        if (lower(metodo)=='e')
            title('Ejercicio 5 - Explicito');
        elseif (lower(metodo) == 'i')
            title('Ejercicio 5 - Implicito');
        elseif (lower(metodo) == 'c')
            title('Ejercicio 5 - Crank-Nicholson');
        end
        
        F(t)=getframe(gcf);
    end
    
    disp('*** FIN ***')
end