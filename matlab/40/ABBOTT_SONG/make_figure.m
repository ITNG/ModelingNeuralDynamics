subplot(111);

tau_plus=10;
tau_minus=10;
K_plus=0.1;
K_minus=2/3*K_plus;

z=[1:1000]/1000*tau_plus*2;
plot(z,K_plus*exp(-z/tau_plus),'-k','Linewidth',8);
hold on;
plot(z,K_plus*exp(-z/tau_plus).*(1-exp(-abs(z)*5/tau_plus)), ...
    '-r','Linewidth',2);
z=-z;
plot(z,-K_minus*exp(z/tau_minus),'-k','Linewidth',8);
plot(z,-K_minus*exp(z/tau_minus).*(1-exp(-abs(z)*5/tau_plus)), ...
    '-r','Linewidth',2);
plot([0,0],[-K_minus*1.2,K_plus*1.2],'--k','Linewidth',1);
axis([-2*tau_plus,2*tau_plus,-K_minus*1.2,K_plus*1.2]);
plot([-2*tau_plus,2*tau_plus],[0,0],'--k','Linewidth',1);
%set(gca,'Xtick',[0]); 
set(gca,'Ytick',[0]);

hold off;
set(gca,'Fontsize',24);
xlabel('$z$ [ms]','Fontsize',32);
title('$F_0$ (black) and $F$ (red)','Fontsize',32);

shg;