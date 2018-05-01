 clc;clf;close all;clear all;
Lx=20;Nx=800; maxerror=1e-10; nmax=2000;k0=2*pi/Lx;
% dx=Lx/Nx; x=-Lx/2:dx:Lx/2-dx; kx=[0:Nx/2-1  -Nx/2:-1]*2*pi/Lx;
dx=Lx/Nx; x=-Lx/2:dx:Lx/2; kx=[0:Nx/2  -Nx/2:-1]*2*pi/Lx;
K2=kx.^2; alpha=6;
g=(sin(x)).^2;
gx=-2*sin(x).*cos(x);
V=g.^2+alpha*g+i*gx;
c=3.8; 
DT=0.7;
%0.5 0.7 0.9 0.1 0.3 1 1.2 1.4 2 1.8 1.6 1.7 4 0.8 1.9 0.6 2
freq=0;
P=1.5;
% for c=10:-0.1:0.1
%     c
%     flag=0;
%  for DT=0.1:0.1:10;%Ä¬ÈÏÖµc=3.8,DT=0.7
%      if flag==1
%          break;
%      end
%     for P=2:0.1:3
        
 u=sech(2*x);                 % initial condition

 u=u*sqrt(P/(sum(abs(u.*conj(u)))*dx));   % power normalization
 for nn=1:nmax                              % iteration starts
    L00u=ifft(-K2.*fft(u))+(u.*conj(u)+V).*u;
    Minvu=ifft(fft(u)./(c+K2));
    mu=sum(Minvu.*L00u)/sum(Minvu.*u);
    L0u=L00u-mu*u;
    uerror(nn)=max(abs(L0u));
    if uerror(nn) < maxerror
        break
    end
    u=u+ifft(fft(L0u)./(c+K2))*DT;
    u=u*sqrt(P/(sum(abs(u.*conj(u)))*dx));
 end                                        % iteration ends
%    if mu>4.5
%     break;
%        flag=1;
%    end
%    if uerror(nn)<maxerror && mu<4.6
% %    sa=[c DT P nn uerror(nn)];
% %      save result.txt sa -ascii -append
%    figure(1)
%    plot(real(mu),P,'b.'),hold on,axis([3.5 4.5 0 3])
%    else
%        flag=1;
%           break;
%    end
%     end
%      end
% end
figure(1)
plot(x, real(V),'k-',x,imag(V),'b--'),axis([-8 8 -3 9])