I=[-100:0]/100;

subplot(111);
plot(I,zeros(101,1),'-k','Linewidth',2);
hold on;

plot(-I,zeros(101,1),'--k','Linewidth',2);
plot(I,sqrt(-I),'--k','Linewidth',2);
plot(I,-sqrt(-I),'--k','Linewidth',2);
set(gca,'Fontsize',20);
xlabel('$I$','Fontsize',26);
ylabel('oscillation amplitude','Fontsize',26)
hold off;
axis([-1,1,-1,1]); axis('square');
shg;