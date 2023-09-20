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
%
% tiempo incial
%
tiempo=0.0;
%
% Valor inicial. 
% u0 array solucion en el paso de tiempo t^n
%
% Talla de u0 es 1x(M+2)
%
% Uso de dato inicial irregular
% Muestra el uso de vectores
%
u0 = zeros(1,M+2);
% for j=2:M+1
%   u0(j)=1;
%   end
for j=1:M+2
    if ((xe(j)>0.3)&&(xe(j)<0.5))
        u0(j)=1;
    end
    if ((xe(j)>0.7)&&(xe(j)<0.9))
        u0(j)=4;
    end
end
%
% Coeficiente de difusion
%
k=0.05;
a_t = -5; % Coeficiente de transporte 
%
% u1 array solucion en el paso de tiempo t^{n+1}
%
u1 = zeros(1,M+2); 
%
% Tiempo final
%
T=0.2;
%
%
Nfin=80;
dt=T/Nfin;
mu=k*dt/h2;
F=a_t*dt/(2.0*h);
%
%    Diagonales de la matriz
%
c = (F-mu)*ones(M,1);  % de 1 a M-1 diagonal superior
b = (1+2.0*mu)*ones(M,1);  % de 1 a M diagonal central
a = (-F-mu)*ones(M,1);  % de 2 a M diagonal inferior
[al,bu]=thomasLUfact(a,b,c);
x=zeros(M,1);
y=zeros(M,1);

%    
% Construccion del sistema MxM
%    Diagonales de la matriz
%
dcent = (1+2.0*mu)*ones(M,1);  % d valores indexados de 1 a M
dinf = (-F-mu)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
dsup = (F-mu)*ones(M-1,1);  % d-1 valores indexados de 1 a M-1
%
% Construimos matriz a partir de las diagonales principales
%
Pim=diag(dinf, -1) + diag(dcent, 0) + diag(dsup, 1);
%
% termino independiente
%
d=zeros(M,1);
%
% Dato inicial
%
plot(xe,u0,'+');
title(['Aprox con M = ',num2str(M),...
       ' Dato inicial ',num2str(tiempo)]);
axis([0 1 0 5]);
%pause;
%
% iteracion temporal
% 
for nt=1:Nfin
 %
 % Avanzamos en nivel de tiempo
 %
 tiempo=nt*dt;
%
%construccion termino independiente
%
% Asignacion de datos al termino 
% independendiente del sistema Ay=d
%
for j=1:M   % rango de 1 a M
    d(j)=u0(j+1);%+dt*f(j+1);
end
    %
%Resolucion con LU
%
% Resolvemos bajada
%
tic
y(1)=d(1);
for i=2:M
y(i)=d(i)-al(i)*y(i-1);
end
%
% Resolvemos subida
%
x(M)=y(M)/bu(M);
for i=M-1:-1:1
x(i)=(y(i)-c(i)*x(i+1))/bu(i);
end
tt=toc;
disp(['Thomas LU  tiempo= ',num2str(tt)]);
%
% Resolucion sistema lineal
%
tic
z = Pim\d; 
tm=toc;
disp(['MATLAB \  tiempo= ',num2str(tm)]);
disp(['MATLAB vs Thomas ',num2str(max(abs(x-z)))]);
disp(['----------------------------- ']);
%
% insertamos la solucion dentro del rango de 1 a M
%
for j=1:M
  u1(j+1)=x(j);
end
% Asignacion de los datos de contorno en x=xIni y en x=xFin
u1(1)=0.0; 
u1(M+2)=0.0; 
%
% Dibujamos la solucion aproximada
%
figure(1);
plot(xe,u1,'.')
title([' kappa = ',num2str(k),...
       ', a = ', num2str(a_t),...
       ', difusion-transporte con dx = ',num2str(h),...
      ' dt ',num2str(dt),...
      ' tiempo ',num2str(tiempo)]);
  axis([0 1 0 5]);
  pause(.1);
%
%actualizamos
%
  for j=1:M+2
   u0(j)=u1(j);
  end
end

  
   

