% Ejercicio 3
clear all;
close all;
a = [0 -5];
h = [0.1 0.05 0.025];
N = [10 20 40];
A = @(a) [a -(1+a) 1];
b = @(a) [-(1+a)/2 (3-a)/2];
k = 2;

f = @(t, y) 4*t.*y.^(1/2);
yt = @(t) (t.^2+1).^2;
t0 = 0;
T = 1;
ht = 0.001;
ttrue = t0:ht:t0+T;
ytrue = yt(ttrue);

for i=1:length(a)
  va = A(a(i));
  vb = b(a(i)); 
  figure(i);
  suptitle_string = sprintf('a=%d', a(i));
  for j=1:length(N)
    t = t0:h(j):t0+(k-1)*h(j);
    y0 = yt(t);
    [t, y] = multipasoExplicito(k, va, vb, f, t0, T, y0, N(j));
    subplot(length(N), 1, j);
    plot(t, y, '*-', ttrue, ytrue, 'r-');
    title_string = sprintf('a=%d, h=%f', a(i), h(j));
    title(title_string);  
    err = max(abs(yt(t)-y));
    printf('a=%d, h=%f, error = %f\n', a(i), h(j), err);
    endfor
endfor

%{
Se puede observar que con a=0 se trata de un método 0-estable y consistente, 
y el error tiende a 0 cuando h tiende a 0, 
mientras que con a = -5, el método es consistente, pero no 0-estable,
y el error se dispara cuando h tiende a 0.
%} 