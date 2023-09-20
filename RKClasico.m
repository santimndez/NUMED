function [t, y] = RKClasico(f, y0, t0, T, N)
  % Parametros del metodo Runge-Kutta clasico
  c = [0 0.5 0.5 1];
  b = [1/6 1/3 1/3 1/6];
  A = [0    0   0   0;
       0.5  0   0   0;
       0    0.5 0   0;
       0    0   1   0];
  [t, y] = RK(c, A, b, f, y0, t0, T, N); %LLama a la funcion general
