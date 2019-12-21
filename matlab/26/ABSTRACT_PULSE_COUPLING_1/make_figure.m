clear; clf;

phi=[0:100]/100;

subplot(221);
plot(phi,g(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$g$','Fontsize',20);
axis([0,1,0,1],'square');

subplot(222);
plot(phi,f(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$f$','Fontsize',20);
axis([0,1,0,1],'square');
f_left=f(phi(1:100));
f_right=f(phi(2:101));
if min(f_right-f_left)<=0,
    disp('f is not strictly increasing');
end;
    

subplot(223);
plot(phi,bigF(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$F$','Fontsize',20);
axis([0,1,0,1],'square');

subplot(224);
plot(phi,bigG(phi),'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$\varphi$','Fontsize',20);
ylabel('$G$','Fontsize',20);
axis([0,1,0,1],'square');
hold on;
plot([0,1],[0,1],'--k','Linewidth',1);
hold off;

shg;