
clear;
clf;

tau_m=10;
I=0.12; 
g_I=0.015; 
tau_I=5;
dt=0.01;
x=(1:9999)'/10000;

subplot(221);
psi_vec=psif(tau_m,I,g_I,tau_I,dt,x);
plot(x,psi_vec,'.k','Markersize',5);
axis([0,1,0,1]);
hold on;
plot([0,1],[0,1],'--k','Linewidth',1);
set(gca,'Fontsize',16);
xlabel('$x$','Fontsize',20);
ylabel('$\psi$','Fontsize',20);
hold off;
axis('square');
shg;

subplot(222);
psi_vec=psif(tau_m,I,g_I,tau_I,dt,x);
phi_vec=psif(tau_m,I,g_I,tau_I,dt,psi_vec);
plot(x,phi_vec,'.k','Markersize',5);
axis([0,1,0,1]);
hold on;
plot([0,1],[0,1],'--k','Linewidth',1);
set(gca,'Fontsize',16);
xlabel('$x$','Fontsize',20);
ylabel('$\phi$','Fontsize',20);
hold off;
axis('square');
shg;
