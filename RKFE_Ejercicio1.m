% RK Fehlberg 4(5)
clear all;
close all;
% Función del problema y'=f(t,y)
f = @(t, y) (0<=t && t<=2).*(-y+t) + (2<t && t<=3).*(-3*y+t);
g = @(t, y) (0<=t && t<2).*(-y+t) + (2<=t && t<=3).*(-3*y+t);
F = @(t, y) -y + t;
G = @(t, y) -3*y + t;
y0 = 1;
t0 = 0;
T1 = 2;
T = 3;
N = 1000000;
% Parámetros del método Runge-Kutta Fehlberg 4(5)
p = 4;
c = [0 0.25 3/8 12/13 1 0.5];
beta = [25/216 0 1408/2565 2197/4104 -1/5 0];
b = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
A = [0          0           0           0         0       0;
     0.25       0           0           0         0       0;
     3/32       9/32        0           0         0       0;
     1932/2197  -7200/2197  7296/2197   0         0       0;
     439/216    -8          3680/513    -845/4104 0       0;
     -8/27      2           -3544/2565  1859/4104 -11/40  0];
tol = 1e-4;
h0 = 0.5;
%% Aplicamos RK Adaptativo a los dos trozos de la funcion
%%[t, yRK] = RKAdaptativo(c, A, beta, b, f, y0, t0, T, h0, N, tol, p);
[t1, yRK1] = RKAdaptativo(c, A, beta, b, f, y0, t0, T1, h0, N, tol, p);
y1 = yRK1(length(t1));
[t2, yRK2] = RKAdaptativo(c, A, beta, b, g, y1, T1, T-T1, h0, N, tol, p);
%% Solución exacta
ttrue = t0:0.001:t0+T; %Particion fina
ytrue = @(t) (0<=t & t<=2).*(t-1+2*exp(-t)) + (2<t & t<=3).*(1/3*t-1/9+(4/9*exp(6)+2*exp(4))*exp(-3*t));
y = ytrue(ttrue);
%% Dibujamos la solución
t = [t1, t2];
y2 = ytrue(t);
yRK = [yRK1, yRK2];
NFin = length(t)-1;
plot(t,yRK,'o-',ttrue,y,'r-');
err=max(abs(yRK-y2));
disp([' Error global  = ',num2str(err)]);
disp([' Nfin  = ',num2str(NFin)]);