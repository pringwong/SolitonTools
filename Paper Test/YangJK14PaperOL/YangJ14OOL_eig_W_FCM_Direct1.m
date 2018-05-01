% p1.m: Fourier collocation method for the eigenvalue spectrum of
% the Zakharov-Shabat system, output shown in Fig. 2.2(a)
clc;clf;close all;clear all;
  N=1500; L=500; Nx=10*N;         % number of Fourier modes is 2N+1
  dx=L/Nx; x=-L/2:dx:L/2; k0=2*pi/L;
  A=2;x0=1.2;alpha=-0.9;
 g=A*(exp(-(x+x0).^2)+exp(-(x-x0).^2));                % potential function
 gx=-2*A*(x+x0).*exp(-(x+x0).^2)-2*A*(x-x0).*exp(-(x-x0).^2);
 V=g.^2+alpha*g+i*gx;
  for n=-N:1:N
     C(n+N+1)=dx*sum(V.*exp(-i*k0*n*x))/L;
  end
  B1=-k0^2*diag([-N:N])^2; 
  B2=toeplitz([C(N+1:2*N+1) zeros(1,N)],[C(N+1:-1:1) zeros(1,N)]); 
  M=B1+B2;
eigvalues=eig(M);          % calculate eigenvalues
% eigv=diag(eigvalues);
% vector
plot(real(eigvalues),imag(eigvalues), 'k.'),axis([-2 5 -0.5 0.5])
eigvalues(find(real(eigvalues)>0))