% CALOR2D
nnodos = size(xnod,1); % cantidad de nodos totales
nelem = size(icone,1); % cantidad de elementos

K = zeros(nnodos,nnodos);
F = zeros(nnodos,1);

% matriz de rigidez
for i = 1:nelem
    f = icone(i,:);  % vector fila (recorro por filas)
    
    % coordenadas de los nodos de un elem generico
    xi = xnod(f(1),1); yi = xnod(f(1),2);
    xj = xnod(f(2),1); yj = xnod(f(2),2);
    xk = xnod(f(3),1); yk = xnod(f(3),2);

    AA = det([1 xi yi; 1 xj yj; 1 xk yk]); % calculo del area

    % calculo de los factores
    ai = (xj*yk - xk*yj)/AA; bi = (yj - yk)/AA; ci = (xk - xj)/AA;
    aj = (xk*yi - xi*yk)/AA; bj = (yk - yi)/AA; cj = (xi - xk)/AA;
    ak = (xi*yj - xj*yi)/AA; bk = (yi - yj)/AA; ck = (xj - xi)/AA;
    
    % matriz de rigidez del elemento lineal 
    KE = [bi^2  + ci^2    bi*bj + ci*cj   bi*bk + ci*ck
          bi*bj + ci*cj   bj^2  + cj^2    bj*bk + cj*ck
          bi*bk + ci*ck   bj*bk + cj*ck   bk^2  + ck^2];
	
	FE = zeros(3,1);
	
	% condiciones Robin (conveccion)
	ncon = size(c_con,1); % cantidad de nodos con condiciones de conveccion
	Kcon = zeros(3,3);
	for j = 1:ncon
		if i==c_con(j,1)
			Kcon(c_con(j,2),c_con(j,2)) = 2;
			Kcon(c_con(j,3),c_con(j,3)) = 2;
			Kcon(c_con(j,2),c_con(j,3)) = 1;
			Kcon(c_con(j,3),c_con(j,2)) = 1;
			Kcon = Kcon*(h*c_con(j,4)/6);
			
			FE(c_con(j,2)) = FE(c_con(j,2)) + h(c_con(j,6))*c_con(j,5)*c_con(j,4)/2;
			FE(c_con(j,3)) = FE(c_con(j,3)) + h(c_con(j,6))*c_con(j,5)*c_con(j,4)/2;
		end
	end
	
	% condiciones Neumann (flujo)
	nneu = size(c_neu,1);
	for j = 1:nneu
		if i==c_neu(j,1)
			FE(c_neu(j,2)) = FE(c_neu(j,2)) - c_neu(j,5)*c_neu(j,4)/2;
			FE(c_neu(j,3)) = FE(c_neu(j,3)) - c_neu(j,5)*c_neu(j,4)/2;
		end
	end
	
	% fuente volumetrica
	nvol = size(f_vol,1);
	for j = 1:nvol
		if i==f_vol(j,1)
			FE = FE + [1; 1; 1]*(G_vol(f_vol(j,2))*t*AA/6);
		end
	end

	% fuente puntual
	npun = size(f_pun,1);
	for j = 1:npun
		if i==f_pun(j,1)
			FE = FE + G_pun(f_pun(j,4))*[ai + bi*f_pun(j,2) + ci*f_pun(j,3);
                                         aj + bj*f_pun(j,2) + cj*f_pun(j,3);
                                         ak + bk*f_pun(j,2) + ck*f_pun(j,3)];
		end
	end
	
    fac = 0.5*k*t*AA; % factor de la matriz K (difusiva)
    K(f,f) = K(f,f) + KE*fac + Kcon; % ensamble
	F(f) = F(f) + FE;
end

% condiciones Dirichlet
ndir = size(c_dir,1);
for j = 1:ndir
	K(c_dir(j,1),:) = 0;
	K(c_dir(j,1),c_dir(j,1)) = 1;
	
	F(c_dir(j,1)) = c_dir(j,2);
end

TEMP = K\F;

% grafica
TempMin = min(min(TEMP));
TempMax = max(max(TEMP));
caxis([TempMin - 1, TempMax + 1]); % pseudocolor axis scaling
xmin = min(xnod(:,1)); 
xmax = max(xnod(:,1));
ymin = min(xnod(:,2)); 
ymax = max(xnod(:,2));
wx = xmax - xmin; 
wy = ymax - ymin;
axes = [xmin - wx*0.1, xmax + wx*0.1, ymin - wy*0.1, ymax + wy*0.1];
axis(axes);
axis('equal');

cla % clear current axis.
fi3 = patch('Faces', icone, 'Vertices', xnod, 'FaceVertexCData', TEMP, 'FaceColor', 'interp', 'EdgeColor', 'k', 'EraseMode', 'normal'); 

% dibuja una marca en los nodos
lin = line(xnod(:,1), xnod(:,2), 'LineStyle', 'none', 'Marker', 'o');
title('Calor2D');
grid on;
colorbar;

% guarda la grafica en un archivo
print('-djpeg', 'calor2D');