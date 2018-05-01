clc;clf;close all;clear all;
Lt=6;Nt=800; maxerror=1e-10; nmax=2000;w0=2*pi/Lt;
% dx=Lx/Nx; x=-Lx/2:dx:Lx/2-dx; kx=[0:Nx/2-1  -Nx/2:-1]*2*pi/Lx;
dt=Lt/Nt; t=-Lt/2:dt:Lt/2; wt=[0:Nt/2  -Nt/2:-1]*w0;

beta2=-0.02;
beta3=0.0030;
gamma=0.21;%0.09
omega=2;%2
alpha=0.0193;%0.0005
g=0.02;%0.0035
TR=0.0045;%0.0045

ln=1;
e=-0;
u=0.9*exp(-t.^2);
% u=0.9*exp(-(t-3).^2).*exp(-i*2*t)+0.9*exp(-(t+3).^2).*exp(i*2*t);
u_p=(1+e)*u;
Sz=8*Nt;
Lz=60;

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
axis([0 Nt 0 Sz])
xlabel('t/fs'); ylabel('z/mm'); zlabel('amplitude');
set(gca, 'xtick', [0 Nt/2 Nt], 'xticklabel',{-Lt/2',0,Lt/2})
set(gca, 'ytick', [0 Sz/2 Sz], 'yticklabel',{0,'',Lz})
hold on;
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
set(gca,'linewidth',1,'Fontsize',20,'Fontname','Times');

figure(2)
plot(op_pulse(floor(((1/3+1/2)/2-0.01)*Sz),:),'Color','k','LineWidth',2);
xlabel('t/fs'); ylabel('|u|'); 
xlim([0,Nt])
ylim([0,1.3])
set(gca, 'xtick', [0 Nt/2 Nt], 'xticklabel',{-Lt/2',0,Lt/2})
hold on;
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
set(gca,'linewidth',1,'Fontsize',20,'Fontname','Times');