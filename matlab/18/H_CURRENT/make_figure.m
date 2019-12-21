v=[-100:50];

subplot(221);
plot(v,r_inf(v),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$v$ [mV]','Fontsize',20); ylabel('$r_\infty(v)$','Fontsize',20);
axis([-100,50,0,1]);

subplot(222);
plot(v,tau_r(v),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$v$ [mV]','Fontsize',20); ylabel('$\tau_r$ [ms]','Fontsize',20);
axis([-100,50,0,1000]);

shg;