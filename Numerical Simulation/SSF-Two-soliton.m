clc; clear all; close all; clf;
h=1000;%采样点数目,步数
to=125e-12;
C=0;
beta_2=-20e-27;%-20e-27
gamma=0.00035736;%0.00035736
mu=0;
Po=0.0064;%
t=-4096e-12:1e-12:4095e-12;
dt=1e-12;
Ao=sqrt(Po/4);
step=0;
Ld=(to^2)/(abs(beta_2));
Ln=1/abs(gamma)/Po;
z=3*Ld;
beta_3=6e-37;%1e-36,6e-37
delta=-3e-14;%小于10-14,7e-14,-6.5*10^(-14),-1.2*10^(-13),-3e-14

for m1=1:1:1
ln=1;
omega=3;
u=Ao*(sech(t/to+step).*exp(-i*omega*t/to)+sech(t/to-step).*exp(i*omega*t/to));%page#47 G.P.AGrawal, not effective enough
l=length(u);
dw=1/l/dt*2*pi;
w=(-1*l/2:1:l/2-1)*dw;
u=fftshift(u);
w=fftshift(w);
spectrum=fft(fftshift(u));
Dw=exp((i*beta_2*w.^2-i*beta_3*w.^3+i*mu)*(h/2));
for jj=0:h:z-mod(z,h)
spectrum=Dw.*spectrum;
f=ifft(spectrum);
nd=delta*gradient(f)./dt./f;%循环边界条件
nd(find(isnan(nd)))=0;
f=f.*exp((i*gamma-nd).*((abs(f)).^2)*(h));
spectrum=fft(f);
spectrum=Dw.*spectrum;
f=ifft(spectrum);
op_pulse(ln,:)=abs(f);%saving output pulse at all intervals
ln=ln+1;
end


n=figure(m1);
mesh(op_pulse(1:1:ln-1,:));
view([0 90]);
title({'Pulse Evolution'});
xlabel('Time'); ylabel('distance'); zlabel('amplitude');
grid on;
hold on;

% set(n,'visible','off');
% saveas(gcf,['omega=',num2str(omega),',b2=20e-27,g=0.00035736,b3=6e-37,delta=-3e-14,mu=0.jpg']);
end