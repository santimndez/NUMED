clear all;
%
% Difusion de dato discontinuo
%
format long;
%
% MATLAB indexa desde 1 los vectores
%
% Dimension de la particion M+1 xe(1),xe(2),...,xe(M+2) donde 
%   xe(1)=xIni, xe(M+2)=xFin y los nodos intermedios son 
%     xe(2),xe(3),...,xe(M+1)
%

M=400;
xIni = 0.0;
xFin =1.0;
h = (xFin - xIni)/(M+1);
%
% Talla de x es 1x(M+2)
%
xe =xIni:h:xFin; % particion talla h, total de M+2 puntos

h2=h*h;
% Valor inicial. 
% u0 = 0;
% u1 = 0;
%
% Coeficiente de difusion
%
k=0.05;
a_t = -1; % Coeficiente de transporte
%f = (1*(0.3 < xe(1:M) & xe(1:M) < 0.5) + 4*(0.7< xe(1:M) & xe(1:M)<0.9))';
f = ones(M, 1);
%
% u array solucion
%
u = zeros(1,M+2); 
%
%    Diagonales de la matriz
%
dsup = (a_t*h-2*k)*ones(M-1,1);  % de 1 a M-1
dcent = (4*k)*ones(M,1);  % de 1 a M
dinf = (-a_t*h-2*k)*ones(M-1,1);  % de 2 a M
d = 2*h2*f; % termino independiente
%
% Construimos matriz a partir de las diagonales principales
%
Pim=diag(dinf, -1) + diag(dcent, 0) + diag(dsup, 1);
z = Pim\d;
for j=1:M
  u(j+1)=z(j);
end 

%
% Dibujamos la solucion aproximada
%
figure(1);
plot(xe,u,'.')
title([' kappa = ',num2str(k),...
       ', difusion con dx = ',num2str(h),...
      'transporte con a = ', num2str(a_t)]);
  axis([0 1 0 5]);
  pause(.1);

  
   

