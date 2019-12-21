clear; clf;
N=200;
phi=[0:N]/N;
tau_m=2; 
I=0.13;

g_hat=1/pi/sqrt(tau_m*I-1/4)./(1+(tan(pi*(phi-1/2))).^2);


subplot(111);
plot(phi,g_hat,'-k','Linewidth',6);
set(gca,'Fontsize',24);
xlabel('$\varphi$','Fontsize',32);
ylabel('$\hat{g}$','Fontsize',32);
axis('square');

shg;