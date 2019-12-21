v=[-100:50]+10^(-3);


subplot(321);
plot(v,m_inf(v),'--k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$m_{\infty}$','Fontsize',18);
axis([-100,50,0,1]);
hold off;

% subplot(322);
% plot(v,tau_m(v),'--k','Linewidth',2);
% set(gca,'Fontsize',14);
% ylabel('$\tau_m$ [ms]','Fontsize',18);
% axis([-100,50,0,0.3]);

subplot(323);
plot(v,h_inf(v),'--k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$h_{\infty}$','Fontsize',18);
axis([-100,50,0,1]);
hold off;

subplot(324);
plot(v,tau_h(v),'--k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$\tau_h$ [ms]','Fontsize',18);
axis([-100,50,0,20]);
hold off;

subplot(325);
plot(v,n_inf(v),'--k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$n_{\infty}$','Fontsize',18);
axis([-100,50,0,1]);
hold off;

subplot(326);
plot(v,tau_n(v),'--k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$\tau_n$ [ms]','Fontsize',18);
axis([-100,50,0,20]);
hold off;

shg;
