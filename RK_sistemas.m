function [t,y] = RK_sistemas(n, c, A, b, f, y0, t0, T, N)
  % Aplica el metodo de Runge-Kutta a y'=f(t,y), con y de dimension n

  % Parametros:
  % c vector suma fila i de A
  % A matriz A del tablero de Butcher para calcular cada pendiente
  % b vector para el metodo que se usa para estimar el error (mismo tamano que beta)

  % f(t,y) funcion que define la EDO, la funcion devuelve un vector columna nx1
  % (t0,y0) condicion inicial, y0 vector columna
  % [t0,t0+T] intervalo de la solucion
  % N numero de puntos de la particion

  % Returns:
  % t vector con los valores de t en cada paso (la particion)
  % y vector con los valores de y en cada paso

    h = T/N;
    t(1) = t0;
    y(:, 1) = y0;
    [m,m] = size(A);
    for i = 1:N
      q = zeros(n, 1); % pendiente promediada para avanzar
      k(:, 1) = f(t(i), y(:, i)); % vector de pendientes
      for j = 2:m
        k(:, j) = f(t(i) + c(j)*h, y(:, i)+h*k(:, 1:j-1)*A(j,1:j-1)'); % Actua igual que el metodo para dimension 1 pero q, y son vectores nx1 y k es una matriz nxm
        q = q + k(:, j)*b(j);
      endfor
      q = q + k(:, 1)*b(1);
      dy = h*q;
      y(:, i+1) =  y(:, i) + dy; % Se avanza
      t(i+1) =  t(i) + h;
    endfor
