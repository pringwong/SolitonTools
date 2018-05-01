% p1.m: Fourier collocation method for the eigenvalue spectrum of
% the Zakharov-Shabat system, output shown in Fig. 2.2(a)
clc;clf;close all;clear all;
  N=500; L=200; Nx=2*N;         % number of Fourier modes is 2N+1
  dx=L/Nx; x=-L/2:dx:L/2-dx; k0=2*pi/L;
  A=2;x0=1.2;alpha=1;
  u=A*(exp(-(x+x0).^2)+exp(-(x-x0).^2))+alpha/2;                % potential function
  for n=-N:1:N
     C(n+N+1)=dx*sum(u.*exp(-i*k0*n*x))/L;
  end
  B1=i*k0*diag([-N:N]); 
  B2=toeplitz([C(N+1:2*N+1) zeros(1,N)],[C(N+1:-1:1) zeros(1,N)]); 
  M=[ -B1   B2
       B2'  B1 ];
  eigvalues=eig(-i*M);          % calculate eigenvalues
  plot(eigvalues, '.')          % plotting spectrum
%   axis([-2 2 -1.7 1.7])
eigvalues(find(abs(imag(eigvalues))>0.6))