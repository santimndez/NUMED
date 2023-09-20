% Ejercicio 2

clear all;
%close all;

f = @(t, w) [w(3); 
             w(4);
             -2*w(1)./(w(1).^2+w(2).^2);
             -2*w(2)./(w(1).^2+w(2).^2)]; 
             % w representa el vector columna [x(t); y(t); v1(t); v2(t)]
t0 = 0;
T = 8;
w0 = [-1; 0; 0.1; -0.1];
N = 4000;
h = T/N;
[t, w] = RKClasico_sistemas(4, f, w0, t0, T, N); % Aplicamos el metodo de Runge-Kutta Clasico para el sistema de dimension 4

figure(1);
plot(w(1, :), w(2, :));
title_string = sprintf('trayectoria (x,y) en [0,8] con h=%f', h);
title(title_string);