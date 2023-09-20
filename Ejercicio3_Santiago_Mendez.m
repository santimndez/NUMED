% Ejercicio 3

clear all;
%%close all;
a = [0 -5];
h = [0.1 0.05 0.025];
N = [40 80 160];
A = @(a) [a -(1+a) 1];
b = @(a) [-(1+a)/2 (3-a)/2];
k = 2;

f = @(t, y) t.^3 -2*t.*y;
yt = @(t) 0.5*(t.^2-1) + exp(1-t.^2);
t0 = 1;
T = 4;
ht = 0.001;
ttrue = t0:ht:t0+T;
ytrue = yt(ttrue);

for i=1:length(a)
  va = A(a(i));
  vb = b(a(i)); 
  figure(i);
  for j=1:length(N)
    t = t0:h(j):t0+(k-1)*h(j);
    y0 = yt(t);
    [t, y] = multipasoExplicito(k, va, vb, f, t0, T, y0, N(j));
    subplot(length(N), 1, j);
    if (i==1)
      plot(t, y, '*-', ttrue, ytrue, 'r-'); % En azul la aproximacion y en rojo la solucion exacta
    else
      semilogy(t, y, '*-', ttrue, ytrue, 'r-'); % En azul la aproximacion y en rojo la solucion exacta
    endif
      legend('Multipaso', 'Solucion exacta', 'Location', 'northwest'); 
    title_string = sprintf('a=%d, h=%f', a(i), h(j));
    title(title_string);  
    err = max(abs(yt(t)-y));
    printf('a=%d, h=%f, error = %f\n', a(i), h(j), err);
    endfor
endfor

%{
Se puede observar que con a=0 se trata de un metodo 0-estable y consistente, 
y el error tiende a 0 cuando h tiende a 0, 
mientras que con a = -5, el metodo es consistente, pero no 0-estable,
y el error se dispara cuando h tiende a 0.
%} 