v=[-100:50]+10^(-3);


subplot(321);
plot(v,alpha_m(v)./(alpha_m(v)+beta_m(v)),'-k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$m_\infty$','Fontsize',18);
axis([-100,50,0,1]);

% subplot(322);
% plot(v,1./(alpha_m(v)+beta_m(v)),'-k','Linewidth',2);
% set(gca,'Fontsize',14); 
% ylabel('$\tau_m$ [ms]','Fontsize',18);
% axis([-100,50,0,1]);

subplot(323);
plot(v,alpha_h(v)./(alpha_h(v)+beta_h(v)),'-k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$h_\infty$','Fontsize',18);
axis([-100,50,0,1]);
	
subplot(324);
plot(v,1./(alpha_h(v)+beta_h(v)),'-k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$\tau_h$ [ms]','Fontsize',18);
axis([-100,50,0,10]);

subplot(325);
plot(v,alpha_n(v)./(alpha_n(v)+beta_n(v)),'-k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$n_\infty$','Fontsize',18);
axis([-100,50,0,1]);
	
subplot(326);
plot(v,1./(alpha_n(v)+beta_n(v)),'-k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$\tau_n$ [ms]','Fontsize',18);
axis([-100,50,0,4]);

shg;
