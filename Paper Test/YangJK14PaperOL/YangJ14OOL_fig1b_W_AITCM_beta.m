      clc;clf;close all;clear all;
Lx=30;Nx=800; maxerror=1e-10; nmax=2000;k0=2*pi/Lx;
% dx=Lx/Nx; x=-Lx/2:dx:Lx/2-dx; kx=[0:Nx/2-1  -Nx/2:-1]*2*pi/Lx;
dx=Lx/Nx; x=-Lx/2:dx:Lx/2; kx=[0:Nx/2  -Nx/2:-1]*2*pi/Lx;
K2=kx.^2;  A=2;x0=1.2;alpha=1;
g=A*(exp(-(x+x0).^2)+exp(-(x-x0).^2));
gx=-2*A*(x+x0).*exp(-(x+x0).^2)-2*A*(x-x0).*exp(-(x-x0).^2);
c=3.8; 
DT=0.5;
%0.5 0.7 0.9 0.1 0.3 1 1.2 1.4 2 1.8 1.6 1.7 4 0.8 1.9 0.6 2
freq=0;
P=0.5;
% for c=10:-0.1:0.1
%     c
%     flag=0;
%  for DT=0.1:0.1:10;%默认值c=3.8,DT=0.7
%      if flag==1
%          break;
%      end
%      for P=0.1:0.1:3
        
%  for beta=0:0.01:1.5
 beta=1.5
V=g.^2+alpha*g+i*beta*gx;
 u=exp(-(x-5).^2);                 % initial condition

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
%    if uerror(nn)<maxerror 
% %        && mu<4.6
% %    sa=[c DT P nn uerror(nn)];
% %      save result.txt sa -ascii -append
%    figure(1)
%    plot(real(mu),P,'g.'),hold on
% %    ,axis([3.5 4.5 0 3])
%    else
%        flag=1;
% %           break;
%    end
%     end
%      end
%  end
 imag(mu)
% [beta imag(mu) imag(mu)/real(mu)]
% if uerror(nn)<maxerror 
% %     figure(1)
% % plot(beta,abs(imag(mu)),'k.');hold on;
% figure(1)
% plot(beta,abs(imag(mu)/real(mu)),'k.');hold on;
% else
figure(1)
plot(beta,mu,'k.');hold on;
% else
%     figure(1)
%     plot(beta,abs(imag(mu)),'g.');hold on;
% end
%  end
% if uerror(nn)<maxerror 
% figure(1)
% plot(x, real(u),'k-',x,imag(u),'b--'),axis([-5 5 -1 1.5])
% end
%   end
% for P=0.1:0.35:3.1; 
%         freq=freq+1;
%  u=sech(2*x);                 % initial condition
%  u=u*sqrt(P/(sum(abs(u.*conj(u)))*dx));   % power normalization
%  for nn=1:nmax                              % iteration starts
%     L00u=ifft(-K2.*fft(u))+(u.*conj(u)+V).*u;
%     Minvu=ifft(fft(u)./(c+K2));
%     mu=sum(Minvu.*L00u)/sum(Minvu.*u);
%     L0u=L00u-mu*u;
%     uerror(nn)=max(abs(L0u));
%     if uerror(nn) < maxerror
%         break
%     end
%     u=u+ifft(fft(L0u)./(c+K2))*DT;
%     u=u*sqrt(P/(sum(abs(u.*conj(u)))*dx));
%   end                                        % iteration ends
% 
% figure(freq)
% plot(x, real(u),'k-',x,imag(u),'b--'),axis([-5 5 -1 1.5])
% title(['P=' num2str(P)])
%  
% end

%FCM测本征值
% SP=-k0^2*diag([-Nx/2:Nx/2])^2;
% Mumm=mu*eye(Nx+1);
%   for n=-Nx/2:1:Nx/2
%      pt(n+Nx/2+1)=dx*sum(V.*exp(-i*k0*n*x))/Lx;
%      ptc(n+Nx/2+1)=dx*sum(conj(V).*exp(-i*k0*n*x))/Lx;
%   end
% u2=u.^2;
% u2conj=conj(u).^2;
% u2mo=u.*conj(u);
%   for n=-Nx/2:1:Nx/2
%      u2f(n+Nx/2+1)=dx*sum(u2.*exp(-i*k0*n*x))/Lx;
%      u2conjf(n+Nx/2+1)=dx*sum(u2conj.*exp(-i*k0*n*x))/Lx;
%      u2mof(n+Nx/2+1)=dx*sum(u2mo.*exp(-i*k0*n*x))/Lx;
%   end
% u2m=toeplitz([u2f(Nx/2+1:Nx+1) zeros(1,Nx/2)],[u2f(Nx/2+1:-1:1) zeros(1,Nx/2)]); 
% u2conjm=toeplitz([u2conjf(Nx/2+1:Nx+1) zeros(1,Nx/2)],[u2conjf(Nx/2+1:-1:1) zeros(1,Nx/2)]); 
% u2mom=toeplitz([u2mof(Nx/2+1:Nx+1) zeros(1,Nx/2)],[u2mof(Nx/2+1:-1:1) zeros(1,Nx/2)]);
% ptm=toeplitz([pt(Nx/2+1:Nx+1) zeros(1,Nx/2)],[pt(Nx/2+1:-1:1) zeros(1,Nx/2)]); 
% ptcm=toeplitz([ptc(Nx/2+1:Nx+1) zeros(1,Nx/2)],[ptc(Nx/2+1:-1:1) zeros(1,Nx/2)]); 
% 
% M=[ -Mumm+SP+ptm+2*u2mom   u2m
%     -u2conjm   Mumm-SP-ptcm-2*u2mom];
% 
% figure(2)
% eigvalues=eig(i*M);
% plot(real(eigvalues),imag(eigvalues), '.');axis([-1,1,-50,50])
% 
% eigvalues(find(eigvalues>0.4))

% figure(2)
% semilogy(1:nn, uerror, 'linewidth', 2);

%传输图
% ln=1;
% e=-1.9;
% u_p=(1+e)*u;
% Sz=3*Nx;
% Lz=3*Lx;
% z=linspace(0,Lz,Sz);
% h=z(2)-z(1);
% Dw=exp((-i*kx.^2)*h/2);
% spectrum=fft(u_p);
% op_pulse(1,:)=abs(u_p);
% for jj=1:Sz-2
%     ln=ln+1;
%     spectrum=Dw.*spectrum;
%     f=ifft(spectrum);
%     f=f.*exp((i*abs(f).^2*h+i*V*h));
%     spectrum=fft(f);
%     spectrum=Dw.*spectrum;
%     f=ifft(spectrum);
%     op_pulse(ln,:)=abs(f);
% end
% figure(3)
% mesh(op_pulse(1:1:ln-1,:));
% view([0,90]);
% title({'Pulse Evolution'});
% axis([0 Nx 0 Sz])
% xlabel('x'); ylabel('z'); zlabel('amplitude');
% set(gca, 'xtick', [0 Nx/2 Nx], 'xticklabel',{-Lx/2',0,Lx/2})
% set(gca, 'ytick', [0 Sz/2 Sz], 'yticklabel',{0,'',Lz})
% grid on;
% hold on;

