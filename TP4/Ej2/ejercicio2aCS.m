clear all
close all
clc

% con script calor2D.m

% coordenadas de los nodos
% xnod=[x y]
xnod=[0 0;
      5 0;
      0 5;
      5 5];

% matriz de conectividades icone
% icone=[i j k]
% filas=elementos
% columnas=numeracion de los nodos
icone=[1 2 3;
       2 4 3];

% condiciones de conveccion
% c_con=[elemento nodo_local nodo_local distancia t_exterior h_pertenece]
% elemento=numero del elemento
% nodo_local=nodos que delimitan la frontera
% distancia=largo de la frontera
% t_exterior=temperatura ambiente
% h_pertenece=elige el valor de h del vector de factores h
c_con=[2 2 3 5 30 1];
h=[1.2];

% condiciones Dirichlet
% c_dir=[nodo temperatura]
% nodo=numero del nodo con temperatura impuesta
% temperatura=valor de la temperatura
c_dir=[2 100;
       4 100];

% condiciones Neumann
% c_neu=[elemento nodo_local nodo_local distancia flujo]
% elemento=numero del elemento
% nodo_local=nodos que delimitan la frontera
% distancia=largo de la frontera
% flujo=flujo de calor
c_neu=[1 1 2 5 0; % aislado
       1 1 3 5 2];

% fuente volumetrica
% f_vol=[elemento g_pertenece]
% elemento=numero del elemento
% g_pertenece=elige el valor de G del vector de fuentes G_vol
f_vol=[2 1];
G_vol=[1.2];

% fuente puntual
% f_pun=[elemento x y g_pertenece]
% elemento=numero del elemento
% x=coordenada en x de la fuente
% y=coordenada en y de la fuente
% g_pertenece=elige el valor de G del vector de fuentes G_pun
f_pun=[1 1 1 1];
G_pun=[5];

% constante del material k
k=2;
% constante del espesor t
t=1;

% llama al script calor2D.m
calor2D;

