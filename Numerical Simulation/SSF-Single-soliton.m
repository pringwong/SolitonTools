clc;clf;close all;clear all;
Lt=400;Nt=800;w0=2*pi/Lt;
% dx=Lx/Nx; x=-Lx/2:dx:Lx/2-dx; kx=[0:Nx/2-1  -Nx/2:-1]*2*pi/Lx;
dt=Lt/Nt; t=-Lt/2:dt:Lt/2; wt=[0:Nt/2  -Nt/2:-1]*w0;

% beta2=0;
% beta3=0;
% gamma=10.58;%0.09
% omega=0.08;%2
% alpha=0;%0.0005，效果明显
% g=0;%0.0035
% TR=0;%0.0045

omega=0.06;
alpha=0.0016;
gamma=0.0018;

w1=-(1-0.6)*omega;

g=9*alpha*omega^2/(9*omega^2-7*w1^2)
TR=-1/(2*w1)
beta2=-6*alpha/(9*omega^2-7*w1^2)
beta3=9*alpha/(-9*omega^2*w1+7*w1^3)

% alpha=0;

ln=1;
e=-0;
u=0.7*exp(-(t/60).^2);
% u=0.9*exp(-(t-3).^2).*exp(-i*2*t)+0.9*exp(-(t+3).^2).*exp(i*2*t);
u_p=(1+e)*u;
Sz=8*Nt;
Lz=1000;

z=linspace(0,Lz,Sz);
h=z(2)-z(1);
Dw=exp(i/2*(beta2+i*g/omega^2)*wt.^2*h/2-i*beta3/6*wt.^3*h/2);
spectrum=fft(u_p);
op_pulse(1,:)=abs(u_p);
for jj=1:Sz-2
    ln=ln+1;
    spectrum=Dw.*spectrum;
    f=ifft(spectrum);
    dne=gradient(abs(f).^2)./dt;%循环边界条件
    dne(find(isnan(dne)))=0;
    f=f.*exp(i*gamma*abs(f).^2*h+(g-alpha)/2*h+i*gamma*TR*dne*h);
    spectrum=fft(f);
    spectrum=Dw.*spectrum;
    f=ifft(spectrum);
    op_pulse(ln,:)=abs(f);
end
figure(1)
mesh(op_pulse(1:1:ln-1,:));
view([0,90]);
axis([0 Nt 0 Sz 0  0.7])
xlabel('t/fs'); ylabel('z/mm');zlabel('|u|');
set(gca, 'xtick', [0 Nt/2 Nt], 'xticklabel',{-Lt/2',0,Lt/2})
set(gca, 'ytick', [0 Sz/2 Sz], 'yticklabel',{0,'',Lz})
set(gca, 'ztick', [0 0.7])
grid off;
hold on;
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
set(get(gca,'ZLabel'),'FontSize',22);
set(gca,'linewidth',1,'Fontsize',20,'Fontname','Times');