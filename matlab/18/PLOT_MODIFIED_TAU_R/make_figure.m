v=[-100:50];

subplot(111);
plot(v,tau_r_modified(v),'--b','Linewidth',2);
hold on;
plot(v,tau_r_original(v),'-k','Linewidth',2);
hold off;
set(gca,'Fontsize',24);
xlabel('$v$ [mV]', 'Fontsize',32);
%title('$\tau_r$ (black) and $q \tau_r$ (blue)','Fontsize',32);
shg;