function [al,bu] = thomasLUfact(a,b,c)
% Factoriza A tridiagonal
% 
% c= diagonal superior  c=[c(1),...,c(n-1),0]
% b= diagonal principal b=[b(1),...,b(n-1),b(n)]
% a= diagonal inferior  a=[0,...,a(n-1),a(n)]
n=length(b);

gam=zeros(n,1);
al=zeros(n,1);
bu=zeros(n,1);

gam(1)=b(1);
for k=2:n
gam(k)=b(k)-a(k)*c(k-1)/gam(k-1);
end
%%
%% Coeficientes L y U
%%
bu(1)=gam(1);
for k=2:n
bu(k)=gam(k);
al(k)=a(k)/gam(k-1);
end


