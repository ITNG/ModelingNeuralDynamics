clear; clf;

T=25; tau=20; epsilon=0.2;

F=@(alpha) (alpha+epsilon)*exp(-T/tau);

subplot(111);
plot([0,1-epsilon],[F(0),F(1-epsilon)],'-k','Linewidth',4);
hold on;
plot([1-epsilon,1],[0,0],'-k','Linewidth',6);
plot([1-epsilon,1-epsilon],[0,F(1-epsilon)],':k','Linewidth',1);
axis([0,1,0,1]);
axis('square');
set(gca,'Fontsize',24);
text(0.5,-0.075,'$\alpha$','Fontsize',32);
ylabel('$F$','Fontsize',32);
set(gca,'Xtick',[0,1]);
plot([0,1],[0,1],':k','Linewidth',2);
hold off;
shg;