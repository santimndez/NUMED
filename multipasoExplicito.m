function [t,y] = multipasoExplicito(k, a, b, f, t0, T, y0, N)
  % Aplica un metodo multipaso a y'=f(t,y)

  % Parametros:
  % k numero de pasos del metodo
  % a vector con los coeficientes a_j (tamano k)
  % b vector con los coeficientes b_j (tamano k-1: al ser explicito, se considera b(k)=0)

  % f(t,y) funcion que define la EDO
  % f debe poder recibir vectores
  % t0 t inicial
  % y0 condicion inicial: es un vector de tamano k
  % [t0,t0+T] intervalo de la solucion
  % N numero de puntos de la particion

  % Returns:
  % t vector con los valores de t en cada paso (la particion)
  % y vector con los valores de y en cada paso
  
  % Trabaja con vectores fila
  
  h = T/N;
  t = t0:h:t0+T;
  y = zeros(1, N+1);
  y(1:k) = y0(1:k); % Inicializa y
  for i=1:N-k+1
    y(i+k) = -a(1:k)*y(i:i+k-1)' + h*b(1:k)*f(t(i:i+k-1), y(i:i+k-1))'; % Aplica la recurrencia
  endfor
end
