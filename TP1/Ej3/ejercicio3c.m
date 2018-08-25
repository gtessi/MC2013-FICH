close all
clear all
clc
% ------------------------------------------------
% datos
% ancho de la placa
Lx=3;

% alto de la placa
Ly=1;

% paso x de la malla
dx=0.5;

% paso y de la malla
dy=0.25;

% condiciones de borde
izq=100;
der=100;
sup=100;
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

% matriz global
K=zeros((Nx+1)*(Ny+1));

% vector de incognitas
%phi=zeros((Nx+1)*(Ny+1),1);

% vector de constantes
b=zeros((Nx+1)*(Ny+1),1);

% ecuacion para la posicion del nodo en la esquina superior derecha
% control
%p=pos((Nx+1),(Ny+1),menor);

% calcula los nodos de la malla
for (i=1:Nx+1)
    for (j=1:Ny+1)
        p=pos(i,j,menor);
        % contornos
        if (i==1)
            % borde izquierdo
            K(p,p)=1;
            b(p)=izq;
        else
            if (i==Nx+1)
                % borde derecho
                K(p,p)=1;
                b(p)=der;
            else
                if (j==1)
                    % borde inferior
                    K(p,p)=-2/(dx*dx)-2/(dy*dy);
                    p2=pos(i+1,j,menor);
                    K(p,p2)=1/dy*dy; % este
                    p2=pos(i-1,j,menor);
                    K(p,p2)=1/dy*dy; % oeste
                    p2=pos(i,j+1,menor);
                    K(p,p2)=1/2*dx*dx; % norte
                    b(p)=2*(x(i)^2+y(j)^2)-qinf*2*dy*dx*dx;
                else
                    if (j==Ny+1)
                        % borde superior
                        K(p,p)=1;
                        b(p)=sup;
                    else
                        K(p,p)=-2/(dx*dx)-2/(dy*dy);
                        p2=pos(i,j+1,menor);
                        K(p,p2)=1/(dy*dy); % norte
                        p2=pos(i,j-1,menor);
                        K(p,p2)=1/(dy*dy); % sur
                        p2=pos(i-1,j,menor);
                        K(p,p2)=1/(dx*dx); % oeste
                        p2=pos(i+1,j,menor);
                        K(p,p2)=1/(dx*dx); % este
                        b(p)=2*(x(i)^2+y(j)^2);
                    end
                end
            end
        end
    end
end

% resuelve el sistema
phi=K\b;

% grafica la solucion aproximada
figure;
phir=reshape(phi,Ny+1,Nx+1);
[X,Y]=meshgrid(x,y);
mesh(X,Y,phir);
surf(X,Y,phir);
xlabel('x');
ylabel('y');
zlabel('T');
% ------------------------------------------------