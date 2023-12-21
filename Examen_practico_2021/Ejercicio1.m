% Ejercicio 1
clear all;
close all;
% Funcion del problema y'=f(t,y)
f = @(t, y) -4*(t.^3).*(y.^2);
y0 = 1/10001;
t0 = -10;
T = 10;

%Parametros del método Runge-Kutta utilizado
c = [0 1/3 2/3];
b = [3/7 1/7 3/7];
A = [0    0   0;
     1/3  0   0;
     1/3  1/3 0];
N0 = 100;
n = 15;

%% Solución exacta
ttrue = t0:0.001:t0+T; %Particion fina
ytrue = @(t) 1./(t.^4+1);
yt = ytrue(ttrue);
dydt = @(t) -4*t.^3./(1+t.^4).^2;

plot(ttrue, yt);

N(1) = N0;
for i = 1:n-1
  N(i+1) = 2*N(i);
endfor

for i=1:n
  [t, yRK] = RK(c, A, b, f, y0, t0, T, N(i));
  y = ytrue(t);
  err(i) = max(abs(yRK-y));
  %figure(1);
  %plot(t, yRK, '*-', ttrue, yt, 'r-');
endfor

figure(2);
plot(log(N),log(err),'*',log(N),-log(N),'-', log(N), -2*log(N), '+',log(N),-3*log(N),'-.');
legend('log(err) Runge-Kutta','-log n','-2log n', '-3log n', 'Location','Best');
title(['Santiago Méndez García: Ejercicio 1, orden 1']);
xlabel('N (número de intervalos de la partición)')