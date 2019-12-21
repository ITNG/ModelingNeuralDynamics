theta=[0:100]/100*2*pi;
x=cos(theta);
y=sin(theta);

subplot(131);
plot(x,y,'-k','Linewidth',2);
axis([-1.5,1.5,-1.5,1.5]); axis('square');
axis off;

hold on;
theta0=-0.4*pi;
x0=cos(theta0);
y0=sin(theta0);
eps=0.15;
fill(x0+eps*x,y0+eps*y,'k');

y0=-y0;
fill(x0+eps*x,y0+eps*y,'w','Linewidth',1);

epsilon=0.1;
theta0=0.4;
x0=cos(theta0); y0=sin(theta0);
v=[sin(theta0+0.1),-cos(theta0+0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=-0.4;
x0=cos(theta0); y0=sin(theta0);
v=[sin(theta0+0.1),-cos(theta0+0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=-2/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=2/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

text(-1.0,1.6,'$I < 1/(4\tau_m)$','Fontsize',16);
hold off;

theta=[0:100]/100*2*pi;
x=cos(theta);
y=sin(theta);
subplot(132);
plot(x,y,'-k','Linewidth',2);
axis([-1.5,1.5,-1.5,1.5]); axis('square');
axis off;

hold on;
theta=[0:100]/100*pi;
fill([1-eps,1+eps,1+eps*cos(theta)],[0,0,-eps*sin(theta)],'k');
fill([1-eps,1+eps,1+eps*cos(theta)],[0,0,+eps*sin(theta)],'w','Linewidth',1);

theta0=-2/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=1/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=-1/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=2/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);


theta0=pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

text(-1,1.6,'$I = 1/(4\tau_m)$','Fontsize',16);
hold off;


theta=[0:100]/100*2*pi;
x=cos(theta);
y=sin(theta);
subplot(133);
plot(x,y,'-k','Linewidth',2);
axis([-1.5,1.5,-1.5,1.5]); axis('square');
axis off;

hold on;
theta0=-2/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=1/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=-1/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=2/3*pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);


theta0=pi;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

theta0=0;
x0=cos(theta0); y0=sin(theta0);
v=-[sin(theta0-0.1),-cos(theta0-0.1)]; width=2; col='-k';
arrow(-1.5,1.5,-1.5,1.5,x0,y0,v,epsilon,width,col);

text(-1.1,1.6,'$I > 1/(4\tau_m)$','Fontsize',16);
hold off;


shg;