clear; clf;
v=[-100:50];
w_inf=1./(1+exp(-(v+35)/10));

tau_w=400./(3.3*exp((v+35)/20)+exp(-(v+35)/20));

subplot(221);
plot(v,w_inf,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$v$ [mV]','Fontsize',20); ylabel('$w_\infty$','Fontsize',20);

subplot(222);
plot(v,tau_w,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$v$ [mV]','Fontsize',20); ylabel('$\tau_w$','Fontsize',20);

shg;