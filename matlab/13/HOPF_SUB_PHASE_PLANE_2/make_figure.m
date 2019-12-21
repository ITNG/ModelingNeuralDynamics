clear; clf;
subplot(331);
plot(0,0,'.k','Markersize',25);
set(gca,'Fontsize',12);
xlabel('$x$','Fontsize',12);
ylabel('$y$','Fontsize',12);
A=-.5; B=.5; C=-.5; D=.5;
axis([A,B,C,D]);
axis('square');

hold on;
I=-0.3;
t_final=5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=0.2;
r(1)=r0;
theta(1)=0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=theta(k)+dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=0.0;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');

r(1)=r0;
theta(1)=0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)-dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)-dt*r_inc;
    theta(k+1)=theta(k)-dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
hold off;
title('$I<-1/4$','Fontsize',10);





clear; 
subplot(332);
I=-0.2;
r0=sqrt(1/2-sqrt(1/4+I));
R0=sqrt(1/2+sqrt(1/4+I));
plot(0,0,'.k','Markersize',25);
set(gca,'Fontsize',12);
xlabel('$x$','Fontsize',12);
ylabel('$y$','Fontsize',12);
A=-1.1; B=1.1; C=-1.1; D=1.1;
axis([A,B,C,D]);
axis('square');
hold on;
th=[0:200]/200*2*pi;
xx=cos(th); yy=sin(th);
plot(r0*xx,r0*yy,'--k','Linewidth',2);
plot(R0*xx,R0*yy,'-k','Linewidth',3);

t_final=5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=0.4;
theta0=0;
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=theta(k)+dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=0.2;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)-dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)-dt*r_inc;
    theta(k+1)=theta(k)-dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);

t_final=5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=0.7;
theta0=pi/2;
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=theta(k)+dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=0.2;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)-dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)-dt*r_inc;
    theta(k+1)=theta(k)-dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);

t_final=5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=1.1;
theta0=5*pi/4;
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=theta(k)+dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=0.;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)-dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)-dt*r_inc;
    theta(k+1)=theta(k)-dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);

set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
hold off;
title('$-1/4<I<0$','Fontsize',10);
shg;






clear; 
subplot(333);
I=0.1;
R0=sqrt(1/2+sqrt(1/4+I));
plot(0,0,'ok','Linewidth',1);
set(gca,'Fontsize',12);
xlabel('$x$','Fontsize',12);
ylabel('$y$','Fontsize',12);
A=-1.3; B=1.3; C=-1.3; D=1.3;
axis([A,B,C,D]);
axis('square');
hold on;
th=[0:200]/200*2*pi;
xx=cos(th); yy=sin(th);
plot(R0*xx,R0*yy,'-k','Linewidth',3);

t_final=5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=1.3;
theta0=3*pi/4;
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=theta(k)+dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=-0.1;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)-dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)-dt*r_inc;
    theta(k+1)=theta(k)-dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);

t_final=5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=0.5;
theta0=pi;
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=theta(k)+dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=0;;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');
r(1)=r0;
theta(1)=theta0;
for k=1:m_steps,
    r_inc=I*r(k)+r(k)^3-r(k)^5;
    r_tmp=r(k)-dt05*r_inc;
    r_inc=I*r_tmp+r_tmp^3-r_tmp^5;
    r(k+1)=r(k)-dt*r_inc;
    theta(k+1)=theta(k)-dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);


set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
hold off;
title('$I>0$','Fontsize',10);
shg;


