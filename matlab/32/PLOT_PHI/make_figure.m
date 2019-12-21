
clear;
clf;

tau_m=10;
I=0.12;
g_I=0.05;
tau_I=5;
dt=0.01;
x=(1:9999)'/10000;

subplot(111);
psi_vec=psif(tau_m,I,g_I,tau_I,dt,x);
phi_vec=psif(tau_m,I,g_I,tau_I,dt,psi_vec);
plot(x,phi_vec,'-k','Linewidth',4);
axis([0,1,0,1]);
hold on;
plot([0,1],[0,1],'--k','Linewidth',2);
set(gca,'Fontsize',24);
xlabel('$x$','Fontsize',32);
ylabel('$\phi(x)$','Fontsize',32);
hold off;
axis('square');
shg;
