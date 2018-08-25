close all
clear all
clc
% ------------------------------------------------
% datos
% largo de la barra
L=1;

% cantidad de nodos
M=2;

% colocacion puntual
[K_cp,f_cp,a_cp,error_cp]=ej1_cp(M,L);

% Galerkin
[K_G,f_G,a_G,error_G]=ej1_galerkin(M,L);
% ------------------------------------------------