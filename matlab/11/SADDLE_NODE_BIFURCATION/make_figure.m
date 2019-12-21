figure('color','w');
subplot(131);
a=0.45;
b=1;

x_min=0.2;
x_max=x_min+3;
y_min=-0.5;
y_max=2.5;
eps=0.125;

t_final=200;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);
x=zeros(m_steps+1,1); y=x;

x_plus=1/(2*a*b)+sqrt(1/(4*a^2*b^2)-1);
y_plus=a*x_plus;
theta=[0:100]/100*2*pi; rad=0.1;
x=x_plus+rad*cos(theta); y=y_plus+rad*sin(theta);
fill(x,y,'k');
hold on;

for ijk=0:10,
x(1)=x_min; y(1)=y_min+ijk*(y_max-y_min)/10;
if ijk==4,
    y(1)=y(1)+0.05;
end;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==4,
        y0=0.5;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
        y0=0.08;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==7,
        y0=1.0;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;

plot(x,y,'-k','Linewidth',1);
hold on;
end;

for ijk=0:10,
x(1)=x_min+ijk*(x_max-x_min)/10; y(1)=y_min;
if ijk==4, 
    x(1)=x(1)-0.05;
end;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==4,
        y0=0;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==4,
        y0=0.58;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==7,
        y0=0.0;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;


for ijk=0:3,
x(1)=x_max; y(1)=y_min+ijk*(y_max-y_min)/10;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==1,
        y0=0.65;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;

for ijk=0:3,
x(1)=x_min+ijk*(x_max-x_min)/10; y(1)=y_max;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==0,
        y0=1.0;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;


plot(x,y,'-k','Linewidth',1);
hold on;
end;


axis([0.2,2.2,-0.5,1.5]); axis('square');
plot([0.2,2.2,2.2,0.2,0.2],[-0.5,-0.5,1.5,1.5,-0.5],'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);

x_minus=1/(2*a*b)-sqrt(1/(4*a^2*b^2)-1);
y_minus=a*x_minus;
x=x_minus+rad*cos(theta); y=y_minus+rad*sin(theta);
fill(x,y,'w');

hold off;
xlabel('$x$','Fontsize',20);
ylabel('$y$','Fontsize',20);
shg;


subplot(132);
a=0.5;
b=1;

x_min=0.2;
x_max=x_min+3;
y_min=-0.5;
y_max=2.5;


dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);
x=zeros(m_steps+1,1); y=x;

x_plus=1/(2*a*b)+sqrt(1/(4*a^2*b^2)-1);
y_plus=a*x_plus;
theta=[0:100]/100*2*pi; 
x=x_plus+rad*cos(theta); y=y_plus+rad*sin(theta);
fill(x,y,'k');
hold on;

for ijk=0:10,
x(1)=x_min; y(1)=y_min+ijk*(y_max-y_min)/10;
if ijk==4,
    y(1)=y(1)+0.05;
end;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==6,
        y0=0.85;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==5,
        y0=0.3;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;

for ijk=0:10,
x(1)=x_min+ijk*(x_max-x_min)/10; y(1)=y_min;
if ijk==4, 
    x(1)=x(1)-0.05;
end;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==6,
        y0=0;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==10,
        y0=0.35;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;


for ijk=0:3,
x(1)=x_max; y(1)=y_min+ijk*(y_max-y_min)/10;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;

for ijk=0:3,
x(1)=x_min+ijk*(x_max-x_min)/10; y(1)=y_max;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==0,
        y0=1;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==0,
        y0=.61;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;


plot(x,y,'-k','Linewidth',1);
hold on;
end;


axis([0.2,2.2,-0.5,1.5]); axis('square');
plot([0.2,2.2,2.2,0.2,0.2],[-0.5,-0.5,1.5,1.5,-0.5],'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);

x_minus=1/(2*a*b)-sqrt(1/(4*a^2*b^2)-1);
y_minus=a*x_minus;
theta=[0:100]/100*pi+0.73*pi;
x=x_minus+rad*cos(theta); y=y_minus+rad*sin(theta);
fill(x,y,'w');
xlabel('$x$','Fontsize',20);
hold off;
shg;

subplot(133);
a=0.55;
b=1;

x_min=0.2;
x_max=x_min+3;
y_min=-0.5;
y_max=2.5;

t_final=200;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);
x=zeros(m_steps+1,1); y=x;


for ijk=0:10,
x(1)=x_min; y(1)=y_min+ijk*(y_max-y_min)/10;
if ijk==4,
    y(1)=y(1)+0.05;
end;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==8,
        y0=1;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            hold on;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==6,
        y0=0.75;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            hold on;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==7,
        y0=0.375;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            hold on;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;

for ijk=0:10,
x(1)=x_min+ijk*(x_max-x_min)/10; y(1)=y_min;
if ijk==4, 
    x(1)=x(1)-0.05;
end;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
    if ijk==5,
        y0=0;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
    if ijk==8,
        y0=0.;
        if (y(k+1)-y0)*(y(k)-y0)<0,
            u=[x(k+1)-x(k);y(k+1)-y(k)];
            u=u/norm(u);
            u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
            u_left=[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_right(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_right(2)],'-k','Linewidth',1);
            plot([x(k+1),x(k+1)-eps*sqrt(5/4)*u_left(1)], ...
                [y(k+1),y(k+1)-eps*sqrt(5/4)*u_left(2)],'-k','Linewidth',1);
        end;
    end;
end;


plot(x,y,'-k','Linewidth',1);
hold on;

end;


for ijk=0:3,
x(1)=x_max; y(1)=y_min+ijk*(y_max-y_min)/10;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
end;

plot(x,y,'-k','Linewidth',1);
hold on;

end;

for ijk=0:3,
x(1)=x_min+ijk*(x_max-x_min)/10; y(1)=y_max;

for k=1:m_steps,
    x_inc=-a*x(k)+y(k);
    y_inc=x(k)^2/(1+x(k)^2)-b*y(k);
    x_tmp=x(k)+dt05*x_inc;
    y_tmp=y(k)+dt05*y_inc;
    x_inc=-a*x_tmp+y_tmp;
    y_inc=x_tmp^2/(1+x_tmp^2)-b*y_tmp;
    x(k+1)=x(k)+dt*x_inc;
    y(k+1)=y(k)+dt*y_inc;
end;


plot(x,y,'-k','Linewidth',1);

end;

axis([0.2,2.2,-0.5,1.5]); axis('square');
plot([0.2,2.2,2.2,0.2,0.2],[-0.5,-0.5,1.5,1.5,-0.5],'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);


xlabel('$x$','Fontsize',20);
hold off;


shg;
    
