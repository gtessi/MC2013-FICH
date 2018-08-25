function [phi,Nx,Ny,X,Y] = FE(Lx,Ly,Lt,dx,dy,dt)
    % ------------------------------------------------
    % forward euler (explicito)
    % ------------------------------------------------
    % condiciones de borde
    izq=10;
    der=10;
    sup=10;
    % flujo inferior nulo
    qinf=0;
    
    % discretiza el dominio
    x=0:dx:Lx;
    y=0:dy:Ly;
    
    % espaciamientos
    Nx=floor(Lx/dx);
    Ny=floor(Ly/dy);
    
    % direccion con menor cantidad de nodos
    if Nx>Ny
        menor=Ny;
    else
        menor=Nx;
    end
    
    % intervalo de tiempo
    t=0:dt:Lt;

    % cantidad de muestras de tiempo
    nt=length(t);

    % matriz con los valores de la malla en cada instante de tiempo
    phi=zeros((Nx+1)*(Ny+1),nt);

    for i=1:Nx+1
        for j=1:Ny+1
            p=pos(i,j,menor);
            % contornos
            if (i==1)
                % borde izquierdo
                phi(p,1)=izq;
            else
                if (i==Nx+1)
                    % borde derecho
                    phi(p,1)=der;
                else
                    if (j==Ny+1)
                        % borde superior
                        phi(p,1)=sup;
                    end
                end
            end
        end
    end
    
    k=2;
    no_est=1;
    tol=0.0001;

    % calcula la variacion de calor a lo largo del tiempo
    while (k<=nt && no_est==1)
        for i=1:Nx+1
            for j=1:Ny+1
                p=pos(i,j,menor);
                % contornos
                if (i==1)
                    % borde izquierdo
                    phi(p,k)=izq;
                else
                    if (i==Nx+1)
                        % borde derecho
                        phi(p,k)=der;
                    else
                        if (j==1)
                            % borde inferior
                            p_izq=pos(i-1,j,menor); % oeste
                            p_der=pos(i+1,j,menor); % este
                            p_arr=pos(i,j+1,menor); % norte
                            phi(p,k)=dt*(100*(x(i))-2*qinf/dy+(phi(p_der,k-1)+phi(p_izq,k-1))/(dx*dx)-2*phi(p,k-1)*(1/(dx*dx)+1/(dy*dy))+2*phi(p_arr,k-1)/(dy*dy))+phi(p,k-1);
                        else
                            if (j==Ny+1)
                                % borde superior
                                phi(p,k)=sup;
                            else
                                p_izq=pos(i-1,j,menor); % oeste
                                p_der=pos(i+1,j,menor); % este
                                p_arr=pos(i,j+1,menor); % norte
                                p_aba=pos(i,j-1,menor); % sur
                                phi(p,k)=dt*(100*(x(i)+y(j))+(phi(p_der,k-1)+phi(p_izq,k-1))/(dx*dx)-2*phi(p,k-1)*(1/(dx*dx)+1/(dy*dy))+(phi(p_arr,k-1)+phi(p_aba,k-1))/(dy*dy))+phi(p,k-1);
                            end
                        end
                    end
                end
            end
        end
        
        % criterio de corte
        if (norm(phi(:,k)-phi(:,k-1),2)/norm(phi(:,k),2)<tol)
            no_est=0;
        end
        
        k=k+1;
    end
    
    [X,Y]=meshgrid(x,y);
end