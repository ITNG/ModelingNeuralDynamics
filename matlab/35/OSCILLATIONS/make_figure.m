clear; clf;
t=[-250:250]/100;

alpha=0.00001;
phi=exp(alpha*sin(pi*t).^2)-1;
phi=phi/mean(phi(2:length(t)));

subplot(311);
plot(t,phi,'-k','Linewidth',2);
alpha=0;
alpha_str=num2str(alpha);
set(gca,'Fontsize',14);
title(['$\alpha \rightarrow~\!$',alpha_str],'Fontsize',18);
axis([-2.5,2.5,0,7]);

alpha=5;
phi=exp(alpha*sin(pi*t).^2)-1;
phi=phi/mean(phi(2:length(t)));

subplot(312);
plot(t,phi,'-k','Linewidth',2);
alpha_str=num2str(alpha);
set(gca,'Fontsize',14);
title(['$\alpha=~\!$',alpha_str],'Fontsize',18);
axis([-2.5,2.5,0,7]);

alpha=10;
phi=exp(alpha*sin(pi*t).^2)-1;
phi=phi/mean(phi(2:length(t)));

subplot(313);
plot(t,phi,'-k','Linewidth',2);
alpha_str=num2str(alpha);
set(gca,'Fontsize',14);
title(['$\alpha=~\!$',alpha_str],'Fontsize',18);
axis([-2.5,2.5,0,7]);


shg;