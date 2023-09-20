function [t,y] = RK(c, A, b, f, y0, t0, T, N)
  % Aplica el metodo de Runge-Kutta a y'=f(t,y)

  % Parametros:
  % c vector suma fila i de A
  % A matriz A del tablero de Butcher para calcular cada pendiente
  % b vector para el metodo que se usa para estimar el error (mismo tamano que beta)

  % f(t,y) funcion que define la EDO
  % (t0,y0) condicion inicial
  % [t0,t0+T] intervalo de la solucion
  % N numero de puntos de la particion

  % Returns:
  % t vector con los valores de t en cada paso (la particion)
  % y vector con los valores de y en cada paso

    h = T/N;
    t(1) = t0;
    y(1) = y0;
    [m,m] = size(A);
    for n = 1:N
      q = 0; % pendiente promediada para avanzar
      k(1) = f(t(n), y(n)); % vector de pendientes
      for j = 2:m
        k(j) = f(t(n) + c(j)*h, y(n)+h*sum(A(j,1:j-1).*k(1:j-1)));
        q = q + k(j)*b(j);
      endfor
      q = q + k(1)*b(1);
      dy = h*q;
      y(n+1) =  y(n) + dy; % Se avanza
      t(n+1) =  t(n) + h;
    endfor
