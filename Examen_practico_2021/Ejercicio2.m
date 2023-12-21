% Ejercicio 2
clear all;
close all;

alpha = pi/4;
semiperiodo = pi*(1+1/4*sin(alpha/2)^2+9/64*sin(alpha/2)^4);

f = @(t, y) [y(2); -sin(y(1))];
t0 = 0;
y0 = [alpha;0];
T = semiperiodo;
N(1)=10;
n = 10;
for i=1:n-1
  N(i+1) = 2*N(i);
endfor

[t, y] = RKClasico_sistemas(2, f, y0, t0, T, N(n));

figure(1);
plot(y(1, :), y(2, :));
title('Santiago Méndez García: Ejercicio 2, plano de fases');
xlabel('y');
ylabel('dy/dt');

figure(2);
plot(t,y(1, :),t,y(2, :), semiperiodo, 0, '*');
legend('y','dy', '(semiperiodo, 0)');



err(n) = abs(y(2, N(n)+1));
for i=1:n-1
  [t, y] = RKClasico_sistemas(2, f, y0, t0, T, N(i));
  err(i) = abs(y(2, N(i)+1));
endfor

for i=1:n
  printf('Error con N=%d: %f\n', N(i), err(i));
endfor