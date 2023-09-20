function [t,y] = RKClasico_sistemas(n, f, y0, t0, T, N)
  % Parametros del método Runge-Kutta clasico
  c = [0 0.5 0.5 1];
  b = [1/6 1/3 1/3 1/6];
  A = [0    0   0   0;
       0.5  0   0   0;
       0    0.5 0   0;
       0    0   1   0];
  [t, y] = RK_sistemas(n, c, A, b, f, y0, t0, T, N); % Llama a la funcion general
