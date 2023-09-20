% Ejercicio 1

clear all;
%%close all;
% Funcion del problema y'=f(t,y)
f = @(t, y) t.^3-2*t.*y;
y0 = 1;
t0 = 1;
T = 4;

%Parametros del metodo Runge-Kutta utilizado
c = [0 1/2 3/4];
b = [2/9 1/3 4/9];
A = [0    0   0;
     1/2  0   0;
     0    3/4 0];
N0 = 50;
n = 10;

%% Solucion exacta
ttrue = t0:0.001:t0+T; %Particion fina
ytrue = @(t) 0.5*(t.^2-1)+exp(1-t.^2);
yt = ytrue(ttrue);

plot(ttrue, yt);

N(1) = N0;
for i = 1:n-1
  N(i+1) = 2*N(i);
endfor

for i=1:n
  [t, yRK] = RK(c, A, b, f, y0, t0, T, N(i)); % Aplicamos el metodo de Runge-Kutta
  y = ytrue(t);
  err(i) = max(abs(yRK-y));
  printf('N = %d, error = %f\n', N(i), err(i));
  if (i>1)
     cociente(i) = err(i-1)/err(i); %% El cociente se aproximara a 2^p, con p el orden del metodo
     printf('error(N=%d)/error(N=%d)=%f\n', N(i-1), N(i), cociente(i));
  endif
  %%figure(2); %% Descomentar para visualizar la aproximacion del metodo
  %%plot(t, yRK, '*-', ttrue, yt, 'r-');
endfor

figure(1);
plot(log(N),log(err),'*',log(N),-log(N),'-', log(N), -2*log(N),log(N),-3*log(N));
legend('log(error con Runge-Kutta)','-log N','-2log N', '-3log N', 'Location','Best');
title(['Ejercicio 1, orden 3']);
xlabel('log(N)')
