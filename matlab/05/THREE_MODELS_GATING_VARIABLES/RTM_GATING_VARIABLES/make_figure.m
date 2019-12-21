clear;
v=[-100:50]+10^(-3);


subplot(321);
plot(v,alpha_m(v)./(alpha_m(v)+beta_m(v)),'-r','Linewidth',2);
set(gca,'Fontsize',14);
%xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$m_\infty(v)$','Fontsize',18);
axis([-100,50,0,1]);
hold on;

% subplot(322);
% plot(v,1./(alpha_m(v)+beta_m(v)),'-r','Linewidth',2);
% set(gca,'Fontsize',14);
% %xlabel('$v$ [mV]','Fontsize',18); 
% ylabel('$\tau_m$ [ms]','Fontsize',18);
% axis([-100,50,0,0.2]);

subplot(323);
plot(v,alpha_h(v)./(alpha_h(v)+beta_h(v)),'-r','Linewidth',2);
set(gca,'Fontsize',14);
%xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$h_\infty(v)$','Fontsize',18);
axis([-100,50,0,1]);
hold on;
	
subplot(324);
plot(v,1./(alpha_h(v)+beta_h(v)),'-r','Linewidth',2);
set(gca,'Fontsize',14);
%xlabel('$v$ [mV]','Fontsize',18); 
ylabel('$\tau_h$ [ms]','Fontsize',18);
axis([-100,50,0,7]);
hold on;

subplot(325);
plot(v,alpha_n(v)./(alpha_n(v)+beta_n(v)),'-r','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); ylabel('$n_\infty(v)$','Fontsize',18);
axis([-100,50,0,1]);
hold on;
	
subplot(326);
plot(v,1./(alpha_n(v)+beta_n(v)),'-r','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$v$ [mV]','Fontsize',18); ylabel('$\tau_n$ [ms]','Fontsize',18);
axis([-100,50,0,2]);
hold on;


shg;
