clear; clf;

v=(-100:50);

subplot(321);
plot(v,a_inf(v),'-k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$a_\infty(v)$','Fontsize',18);
axis([-100,50,0,1]);

subplot(322);
plot(v,tau_a(v),'-k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$\tau_a$ [ms]','Fontsize',18);
axis([-100,50,0,10]);

subplot(323);
plot(v,b_inf(v),'-k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); ylabel('$b_\infty(v)$','Fontsize',18);
axis([-100,50,0,1]);

subplot(324);
plot(v,tau_b(v),'-k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); ylabel('$\tau_b$ [ms]','Fontsize',18);
axis([-100,50,0,650]);


shg;