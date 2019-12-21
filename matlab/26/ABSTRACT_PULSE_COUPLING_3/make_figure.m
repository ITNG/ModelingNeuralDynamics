clear; clf;

phi=[0:100]/100;

subplot(221);
plot(phi,g(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$g$','Fontsize',20);
axis([0,1,-1/2,1/2],'square');



subplot(222);
plot(phi,bigG(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$G$','Fontsize',20);
axis([0,1,0,1],'square');
hold on;
plot([0,1],[0,1],'--k','Linewidth',1);
hold off;

shg;