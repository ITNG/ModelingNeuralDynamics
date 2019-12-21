clear; clf;


phi=(0:1000)/1000;

subplot(221);
plot(phi,g(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$g$','Fontsize',20);
axis([0,1,-1,0],'square');



subplot(222);
plot(phi,bigG(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$G$','Fontsize',20);
axis([0,1,0,1],'square');
hold on;
plot([0,1],[0,1],'--k','Linewidth',1);

phi_left=phi(1:1000);
phi_right=phi(2:1001);
G_left=bigG(phi_left);
G_right=bigG(phi_right);
ind=find((G_left-phi_left).*(G_right-phi_right)<=0);
plot(phi_left(ind(1)),G_left(ind(1)),'.g','Markersize',30);
plot(phi_right(ind(5)),G_right(ind(5)),'.g','Markersize',30);
plot(phi_left(ind(2)),G_left(ind(2)),'or','Markersize',8,'Linewidth',2);
plot(phi_right(ind(4)),G_right(ind(4)),'or','Markersize',8,'Linewidth',2);
plot((phi_left(ind(3))+phi_right(ind(3)))/2, ...
    (G_left(ind(3))+G_right(ind(3)))/2,'.b','Markersize',30);

hold off;

shg;