% p1.m: Fourier collocation method for the eigenvalue spectrum of
% the Zakharov-Shabat system, output shown in Fig. 2.2(a)
clc;clf;close all;clear all;
  N=200; L=100; Nx=2*N;         % number of Fourier modes is 2N+1
  dx=L/Nx; x=-L/2:dx:L/2-dx; k0=2*pi/L;
  A=2;x0=1.2;alpha=1;
  
  for P=0.1:0.2:3;
%   u=A*(exp(-(x+x0).^2)+exp(-(x-x0).^2))+alpha/2;                % potential function
%   for n=-N:1:N
%      C(n+N+1)=dx*sum(u.*exp(-i*k0*n*x))/L;
%   end
  B1=-k0^2*diag([-N:N])^2; 
  B2=P*diag(ones(2*N+1,1));
  M=B1+B2;
  eigvalues=eig(M)          % calculate eigenvalues
    plot(P,eigvalues, '.')          % plotting spectrum
%       mesh(P,eigvalues)  
%   axis([-2 2 -1.7 1.7])
% eigvalues(find(real(eigvalues)>0))
  end