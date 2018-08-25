function [phi,Nx,Ny,X,Y] = BE(Lx,Ly,Lt,dx,dy,dt)
    % ------------------------------------------------
    % backward euler (implicito)
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
    
    % matriz global
    K=zeros((Nx+1)*(Ny+1));
    
    % vector de constantes
    f=zeros((Nx+1)*(Ny+1),1);

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
                    K(p,p)=1;
                    f(p)=izq;
                else
                    if (i==Nx+1)
                        % borde derecho
                        K(p,p)=1;
                        f(p)=der;
                    else
                        if (j==1)
                            % borde inferior
                            K(p,p)=-2/(dx*dx)-2/(dy*dy)-1/dt;
                            p_izq=pos(i-1,j,menor); % oeste
                            K(p,p_izq)=1/(dx*dx);
                            p_der=pos(i+1,j,menor); % este
                            K(p,p_der)=1/(dx*dx);
                            p_arr=pos(i,j+1,menor); % norte
                            K(p,p_arr)=2/(dy*dy);
                            f(p)=-100*x(i)+2*qinf/dy-phi(p,k-1)/dt;
                        else
                            if (j==Ny+1)
                                % borde superior
                                K(p,p)=1;
                                f(p)=sup;
                            else
                                K(p,p)=-2/(dx*dx)-2/(dy*dy)-1/dt;
                                p_izq=pos(i-1,j,menor); % oeste
                                K(p,p_izq)=1/(dx*dx);
                                p_der=pos(i+1,j,menor); % este
                                K(p,p_der)=1/(dx*dx);
                                p_arr=pos(i,j+1,menor); % norte
                                K(p,p_arr)=1/(dy*dy);
                                p_aba=pos(i,j-1,menor); % sur
                                K(p,p_aba)=1/(dy*dy);
                                f(p)=-100*(x(i)+y(j))-phi(p,k-1)/dt;
                            end
                        end
                    end
                end
            end
        end
        
        phi(:,k)=K\f;
        
        % criterio de corte
        if (norm(phi(:,k)-phi(:,k-1),2)/norm(phi(:,k),2)<tol)
            no_est=0;
        end
        
        k=k+1;
    end
    
    [X,Y]=meshgrid(x,y);
end