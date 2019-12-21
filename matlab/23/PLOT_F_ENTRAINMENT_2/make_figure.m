clear; clf;

T=25; tau=50; epsilon=0.6;

F=@(alpha) (alpha+epsilon)*exp(-T/tau);

subplot(111);
plot([0,1-epsilon],[F(0),F(1-epsilon)],'-k','Linewidth',4);
hold on;
plot([1-epsilon,1],[0,0],'-k','Linewidth',8);
axis([0,1,0,1]);
axis('square');
set(gca,'Fontsize',24);
ylabel('$F$','Fontsize',32);
plot([0,1],[0,1],':k','Linewidth',2);

alpha(1)=0;
alpha(2)=F(alpha(1));
plot([0,0],[0,alpha(2)],'-r','Linewidth',2);
plot([0,alpha(2)],[alpha(2),alpha(2)],'-r','Linewidth',1);
alpha(3)=F(alpha(2));
plot([alpha(2),alpha(2)],[alpha(2),alpha(3)],'-r','Linewidth',1);
plot([alpha(2),alpha(3)],[alpha(3),alpha(3)],'-r','Linewidth',1);
plot([alpha(3),alpha(3)],[alpha(3),0],'-r','Linewidth',1);
plot([0,alpha(3)],[0,0],'-r','Linewidth',2);
plot([alpha(2),alpha(2)],[0,alpha(2)],':r','Linewidth',1);
plot(alpha(1),0,'.r','Markersize',30);
plot(alpha(2),0,'.r','Markersize',30);
plot(alpha(3),0,'.r','Markersize',30);

x=0; y=alpha(2)/2; v=[0; 1]; epsilon=0.03; width=2; col='-r';
arrow(0,1,0,1,x,y,v,epsilon,width,col)
x=alpha(2)/2; y=alpha(2); v=[1;0]; epsilon=0.03; width=2; col='-r';
arrow(0,1,0,1,x,y,v,epsilon,width,col)
x=alpha(2); y=alpha(2)+1/2*(alpha(3)-alpha(2)); v=[0;1]; epsilon=0.03; width=2; col='-r';
arrow(0,1,0,1,x,y,v,epsilon,width,col)
x=(alpha(2)+alpha(3))/2; y=alpha(3); v=[1;0]; epsilon=0.03; width=2; col='-r';
arrow(0,1,0,1,x,y,v,epsilon,width,col)
x=alpha(3); y=alpha(3)/2; v=[0;-1]; epsilon=0.03; width=2; col='-r';
arrow(0,1,0,1,x,y,v,epsilon,width,col)
x=alpha(3)*0.4; y=0; v=[-1;0]; epsilon=0.03; width=2; col='-r';
arrow(0,1,0,1,x,y,v,epsilon,width,col)

text(-0.02,-0.08,'$\alpha_1$','Fontsize',32,'Color','r');
text(alpha(2)-0.02,-0.08,'$\alpha_2$','Fontsize',32,'Color','r');
text(alpha(3)-0.02,-0.08,'$\alpha_3$','Fontsize',32,'Color','r');


set(gca,'Xtick',[1]);
set(gca,'Ytick',[0 0.5 1]);
hold off;
shg;