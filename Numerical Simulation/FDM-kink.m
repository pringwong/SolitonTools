clc;clf;clear all;close all;
Lt=40;Nt=800;dt=Lt/Nt;t=[-Lt/2:dt:Lt/2-dt]';w=[0:Nt/2-1 -Nt/2:-1]'*2*pi/Lt;
xmax=20;dx=0.02;w2=w.^2;nmax=round(xmax/dx);

b1=0.22;
b2=0.097;
alpha=0.39;
pi=0.000003;
qr=-0.00002;
gammai=1;

pr=2*(b2-b1*alpha)*pi/b1
sigma=-qr/(b1^2*(-2*pr+4*alpha^2*pr+8*alpha*pi))
qi=8*b1^2*alpha*sigma*pr+3*b1^2*sigma*pi-4*b1^2*alpha^2*sigma*pi
gammar=-(2*b1*alpha-b2)^2*pi
cr=-sigma*qr
ci=-sigma*qi

u=0.5./(1+exp(t));
  
udata=u; xdata=0;
  for nn=1:nmax                               % integration begins
    du1=-(i*pr-pi)*ifft(-w2.*fft(u))+(gammar+i*gammai)*u+(i*qr-qi)*abs(u).^2.*u+(i*cr-ci)*abs(u).^4.*u;v=u+0.5*du1*dt;
    du2=-(i*pr-pi)*ifft(-w2.*fft(v))+(gammar+i*gammai)*v+(i*qr-qi)*abs(v).^2.*v+(i*cr-ci)*abs(v).^4.*v;v=u+0.5*du2*dt;
    du3=-(i*pr-pi)*ifft(-w2.*fft(v))+(gammar+i*gammai)*v+(i*qr-qi)*abs(v).^2.*v+(i*cr-ci)*abs(v).^4.*v;v=u+du3*dt;
    du4=-(i*pr-pi)*ifft(-w2.*fft(v))+(gammar+i*gammai)*v+(i*qr-qi)*abs(v).^2.*v+(i*cr-ci)*abs(v).^4.*v;
    u=u+(du1+2*du2+2*du3+du4)*dx/6;
     if mod(nn,round(nmax/25)) == 0
        udata=[udata u]; xdata=[xdata nn*dx];
     end
  end                                         % integration ends
wdata=abs(udata');
mesh(t(20:Lt/dt-20), xdata, wdata(:,20:Lt/dt-20));           % solution plotting
% axis([-18 18 0 20])
xlabel('t'); ylabel('x'); zlabel('|u|^2');
% set(gca, 'xtick', [0 Nt/2 Nt], 'xticklabel',{-Lt/2',0,Lt/2})
% set(gca, 'ytick',[300 dex/Lx*Nx Nx/1.4], 'yticklabel',{8,0,-8})
% set(gca, 'ztick', [0 0.05 0.1], 'zticklabel',{0,0.05,0.1})
hold on;
grid off;
view(-37,30);
%   colormap([0 0 0]); view(10, 20)
%   text(-2,  -6, 't', 'fontsize', 15)
%   text(50, 5, 'x', 'fontsize', 15)
%   zlabel('|u|', 'fontsize', 15)
%   axis([-Lt/2 Lt/2 0 xmax]); grid off
%   set(gca, 'xtick', [-40 -20 0 20 40])
%   set(gca, 'ytick', [0 10 20])