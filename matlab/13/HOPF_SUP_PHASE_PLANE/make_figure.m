clear; clf;
subplot(221);
plot(0,0,'.k','Markersize',25);
set(gca,'Fontsize',16);
xlabel('$x$','Fontsize',20);
ylabel('$y$','Fontsize',20);
A=-.5; B=.5; C=-.5; D=.5;
axis([A,B,C,D]);
axis('square');

hold on;
I=-0.02;
t_final=5;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
r(1)=0.2;
theta(1)=0;
for k=1:m_steps,
    r_inc=I*r(k)-r(k)^3;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp-r_tmp^3;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=k*dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
u=[x(2)-x(1); y(2)-y(1)];
phi=0.1;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');


r(1)=0.2;
theta(1)=0;
for k=1:m_steps,
    r_inc=-(I*r(k)-r(k)^3);
    r_tmp=r(k)+dt05*r_inc;
    r_inc=-(I*r_tmp-r_tmp^3);
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=-k*dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
hold off;

title('$I<0$','Fontsize',20);
clear;

subplot(222);
I=0.5;
plot(0,0,'ok','Linewidth',1);
set(gca,'Fontsize',16);
xlabel('$x$','Fontsize',20);
ylabel('$y$','Fontsize',20);
A=-2; B=2; C=-2; D=2;
axis([A,B,C,D]);
axis('square');
hold on;
theta=[0:100]/100*2*pi;
xx=cos(theta); yy=sin(theta);
xx=sqrt(I)*xx; yy=sqrt(I)*yy;
plot(xx,yy,'-k','Linewidth',4);



t_final=3;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=0.4;
r(1)=r0;
theta(1)=0;
for k=1:m_steps,
    r_inc=I*r(k)-r(k)^3;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp-r_tmp^3;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=k*dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
hold on;
u=[x(1000)-x(999); y(1000)-y(999)];
phi=0.3;
R=[cos(phi) sin(phi); -sin(phi) cos(phi)];
u=R*u;
arrow(A,B,C,D,x(1000),y(1000),u,0.075,1,'-k');

t_final=1;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
z=zeros(m_steps+1,1);
r=z; theta=z;

r(1)=r0;
theta(1)=0;
for k=1:m_steps,
    r_inc=-(I*r(k)-r(k)^3);
    r_tmp=r(k)+dt05*r_inc;
    r_inc=-(I*r_tmp-r_tmp^3);
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=-k*dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);

t_final=3;
dt=0.0001; dt05=dt/2;
m_steps=round(t_final/dt);
r0=1.2;
r(1)=r0;
theta(1)=-pi/2;
for k=1:m_steps,
    r_inc=I*r(k)-r(k)^3;
    r_tmp=r(k)+dt05*r_inc;
    r_inc=I*r_tmp-r_tmp^3;
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=-pi/2+k*dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
hold on;
u=[x(2)-x(1); y(2)-y(1)];
arrow(A,B,C,D,x(1),y(1),u,0.075,1,'-k');

t_final=3;
dt=0.0001; dt05=dt/2;
m_steps=round(t_final/dt);
r(1)=r0;
theta(1)=-pi/2;
for k=1:m_steps,
    r_inc=-(I*r(k)-r(k)^3);
    r_tmp=r(k)+dt05*r_inc;
    r_inc=-(I*r_tmp-r_tmp^3);
    r(k+1)=r(k)+dt*r_inc;
    theta(k+1)=-pi/2-k*dt;
end;
x=r.*cos(theta); y=r.*sin(theta);
plot(x,y,'-k','Linewidth',1);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('$I>0$','Fontsize',20);


hold off;





shg;
