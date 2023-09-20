function [t,y] = RKAdaptativo(c, A, beta, b, f, y0, t0, T, h0, N, tol, p)
  % Aplica el m3todo de Runge-Kutta adaptativo a y'=f(t,y)

  % Parametros:
  % c vector suma fila i de A
  % A matriz A del tablero de Butcher para calcular cada pendiente
  % beta vector con pesos para las pendientes para el avance
  % b vector para el método que se usa para estimar el error (mismo tamano que beta)
  % p orden del método de avance

  % f(t,y) función que define la EDO
  % (t0,y0) condición inicial
  % [t0,t0+T] intervalo de la solucion
  % N número de puntos de la particionn

  % Returns:
  % t vector con los valores de t en cada paso (la particion)
  % y vector con los valores de y en cada paso

    h = h0;
    t(1) = t0;
    Tfin = T + t0;
    y(1) = y0;
    n = 1;
    [m,m] = size(A);
    while((n <= N) && (t(n) < Tfin) && (h > eps))
      % RK para calcular pendientes
      q = 0; % pendiente promediada para avanzar
      k(1) = f(t(n), y(n)); % vector de pendientes
      for j = 2 : m
        k(j) = f(t(n) + c(j)*h, y(n)+h*sum(A(j,1:j-1).*k(1:j-1)));
        q = q + k(j)*beta(j);
      endfor
      q = q + k(1)*beta(1);
      dy = h*q;

      tau = abs(sum((b-beta).*k));
      % Aceptamos el paso dado y avanzamos
      if  (tau <= tol)
        y(n+1)=  y(n) + dy;
        t(n+1)=  t(n) + h;

        hNext=min(10*h,Tfin-t(n+1));
        n=n+1;
      % Fin de trabajo cuando se acepta el paso
      else
        %% hNext =h*s
        %% para ser usado si el error local de truncatura no es
        %% lo suficientemente pequeno
        %%
        hNext = 0.9*(tol/tau)^(1/p)*h;
        % Controlamos hNext
     end
    if (hNext< 0.1*h) % Que no sea menor que 0.1*h
      hNext = 0.1*h;
    end
    if (hNext > 10*h) % Que no sea mayor que 10*h
      hNext=10*h;
    end
    if (t(n) + hNext > Tfin)% Que no se sobrepase  t0+T
      hNext = Tfin - t(n);
      Nfin=n;
    end
    h = hNext;
    end

    if ((t(n) < Tfin) && (h <= eps))
      display("Error: no se pudo garantizar la tolerancia con RK");
    endif

end

