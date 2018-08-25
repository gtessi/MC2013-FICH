close all
clear all
clc
% ------------------------------------------------
% datos
% ancho de la placa
Lx=2;

% alto de la placa
Ly=1;

% paso x de la malla
dx=0.5;

% paso y de la malla
dy=0.25;

% duracion
Lt=1;

% delta de tiempo
dt=0.01;

% corre las animaciones para los diferentes metodos
% forward euler (explicito)
movie('e',Lx,Ly,Lt,dx,dy,dt);
% backward euler (implicito)
movie('i',Lx,Ly,Lt,dx,dy,dt);
% crank-nicolson (semi-explicito)
movie('c',Lx,Ly,Lt,dx,dy,dt);
% ------------------------------------------------