clf;

w_max=0.2; tau_w=20; % model parameters
dt=0.0001; dt05=dt/2; % numerical parameters 
epsilon=0.06; width=1; col='k'; % plotting parameters; 

I=-0.15; 

subplot(221);

if I<=0,
    theta_plus=acos((1+I)/(1-I)); theta_minus=-theta_plus;
else
    theta_plus=10; theta_minus=-10;
end;

plot(theta_plus,0,'ok','Markersize',6,'Markerfacecolor','w')
hold on;
plot(theta_minus,0,'ok','Markersize',6,'Markerfacecolor','k')
A=-pi; B=pi; C=0; D=2*w_max;
axis([A,B,C,D],'square');
set(gca,'Fontsize',12);
%xlabel('$\theta$','Fontsize',16);
ylabel('$z$','Fontsize',16);
title('$I=-0.15<I_\ast$','Fontsize',13);
text(-pi-2,0.4,'A','Fontsize',20,'Fontweight','bold');

for ijk=1:25,
    
    theta(1)=-pi; w(1)=w_max*ijk/20*2; k=1;
    a=1.014;
    if ijk>10,
        w(1)=w(1)*a;
    end;

    while theta(k)<pi & w(k)>10^(-6),
        theta_inc=1-cos(theta(k))+(I+w(k))*(1+cos(theta(k)));
        w_inc=-w(k)/tau_w;
        theta_tmp=theta(k)+dt05*theta_inc;
        w_tmp=w(k)+dt05*w_inc;
        theta_inc=1-cos(theta_tmp)+(I+w_tmp)*(1+cos(theta_tmp));
        w_inc=-w_tmp/tau_w;
        theta(k+1)=theta(k)+dt*theta_inc;
        w(k+1)=w(k)+dt*w_inc;
        k=k+1;
    end;
    theta=theta(1:k); w=w(1:k);
    
    plot(theta,w,'-k','Linewidth',1);
    if ijk==10,
        plot(theta,w,'-k','Linewidth',3);
    end;
    
    if ijk==10,
        ind=find(theta>-2);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;
    if ijk==10,
        ind=find(w<0.06);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;
    if ijk==13,
        ind=find(theta>2);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    if ijk==18,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    if ijk==5,
        ind=find(theta>-2);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
end;


x=-2; y=0; v=[1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
x=0; y=0; v=[-1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
x=2; y=0; v=[1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
hold off;

I=-0.1069150434; % very carefully chosen to be close to I_ast

subplot(222);

if I<=0,
    theta_plus=acos((1+I)/(1-I)); theta_minus=-theta_plus;
else
    theta_plus=10; theta_minus=-10;
end;

plot(theta_plus,0,'ok','Markersize',6,'Markerfacecolor','w')
hold on;
plot(theta_minus,0,'ok','Markersize',6,'Markerfacecolor','k')
A=-pi; B=pi; C=0; D=2*w_max;
axis([A,B,C,D],'square');
set(gca,'Fontsize',12);
%xlabel('$\theta$','Fontsize',16);
%ylabel('$z$','Fontsize',16);
title('$I>I_\ast$ but $I \approx I_\ast$','Fontsize',13);
text(-pi-2,0.4,'B','Fontsize',20,'Fontweight','bold');

for ijk=1:25,
    
    theta(1)=-pi; w(1)=w_max*ijk/20*2; k=1;

    while theta(k)<pi & w(k)>10^(-6),
        theta_inc=1-cos(theta(k))+(I+w(k))*(1+cos(theta(k)));
        w_inc=-w(k)/tau_w;
        theta_tmp=theta(k)+dt05*theta_inc;
        w_tmp=w(k)+dt05*w_inc;
        theta_inc=1-cos(theta_tmp)+(I+w_tmp)*(1+cos(theta_tmp));
        w_inc=-w_tmp/tau_w;
        theta(k+1)=theta(k)+dt*theta_inc;
        w(k+1)=w(k)+dt*w_inc;
        k=k+1;
    end;
    theta=theta(1:k); w=w(1:k);
    
    plot(theta,w,'-k','Linewidth',1);
    if ijk==10,
        plot(theta,w,'-k','Linewidth',3);
    end;
    
    if ijk==10,
        ind=find(w<0.05);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;
    
    if ijk==10,
        ind=find(theta>-1.5);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;
    
    if ijk==10,
        ind=find(theta>3);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;
    
    if ijk==5,
        ind=find(theta>-2);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==5,
        ind=find(w<0.02);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==15,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==20,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
end;

x=-2; y=0; v=[1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
x=0; y=0; v=[-1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
x=2; y=0; v=[1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
hold off;

I=-0.05; 
subplot(223);

if I<=0,
    theta_plus=acos((1+I)/(1-I)); theta_minus=-theta_plus;
else
    theta_plus=10; theta_minus=-10;
end;

plot(theta_plus,0,'ok','Markersize',6,'Markerfacecolor','w')
hold on;
plot(theta_minus,0,'ok','Markersize',6,'Markerfacecolor','k')
A=-pi; B=pi; C=0; D=2*w_max;
axis([A,B,C,D],'square');
set(gca,'Fontsize',12);
xlabel('$\theta$','Fontsize',16);
ylabel('$z$','Fontsize',16);
title('$I_\ast < I= -0.05 < I_c=0$','Fontsize',13);
text(-pi-2,0.4,'C','Fontsize',20,'Fontweight','bold');

for ijk=1:25,
    
    theta(1)=-pi; w(1)=w_max*ijk/20*2; k=1;

    while theta(k)<pi & w(k)>10^(-6),
        theta_inc=1-cos(theta(k))+(I+w(k))*(1+cos(theta(k)));
        w_inc=-w(k)/tau_w;
        theta_tmp=theta(k)+dt05*theta_inc;
        w_tmp=w(k)+dt05*w_inc;
        theta_inc=1-cos(theta_tmp)+(I+w_tmp)*(1+cos(theta_tmp));
        w_inc=-w_tmp/tau_w;
        theta(k+1)=theta(k)+dt*theta_inc;
        w(k+1)=w(k)+dt*w_inc;
        k=k+1;
    end;
    theta=theta(1:k); w=w(1:k);
    
    plot(theta,w,'-k','Linewidth',1);
    if ijk==10,
        plot(theta,w,'-k','Linewidth',3);
    end;
    
    if ijk==10,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;

    
    if ijk==5,
        ind=find(theta>-2);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==6,
        ind=find(theta>2);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==15,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==20,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
end;

x=-2; y=0; v=[1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
x=-theta_plus/3; y=0; v=[-1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
x=2; y=0; v=[1; 0]; arrow(A,B,C,D,x,y,v,epsilon,width,col);
hold off;

I=0.05; 
subplot(224);

if I<=0,
    theta_plus=acos((1+I)/(1-I)); theta_minus=-theta_plus;
else
    theta_plus=10; theta_minus=-10;
end;

plot(theta_plus,0,'ok','Markersize',6,'Markerfacecolor','w')
hold on;
plot(theta_minus,0,'ok','Markersize',6,'Markerfacecolor','k')
A=-pi; B=pi; C=0; D=2*w_max;
axis([A,B,C,D],'square');
set(gca,'Fontsize',12);
xlabel('$\theta$','Fontsize',16);
%ylabel('$z$','Fontsize',16);
title('$I=0.05>I_c=0$','Fontsize',13);
text(-pi-2,0.4,'D','Fontsize',20,'Fontweight','bold');

for ijk=1:25,
    
    theta(1)=-pi; w(1)=w_max*ijk/20*2; k=1;

    while theta(k)<pi & w(k)>10^(-6),
        theta_inc=1-cos(theta(k))+(I+w(k))*(1+cos(theta(k)));
        w_inc=-w(k)/tau_w;
        theta_tmp=theta(k)+dt05*theta_inc;
        w_tmp=w(k)+dt05*w_inc;
        theta_inc=1-cos(theta_tmp)+(I+w_tmp)*(1+cos(theta_tmp));
        w_inc=-w_tmp/tau_w;
        theta(k+1)=theta(k)+dt*theta_inc;
        w(k+1)=w(k)+dt*w_inc;
        k=k+1;
    end;
    theta=theta(1:k); w=w(1:k);
    
    plot(theta,w,'-k','Linewidth',1);
    if ijk==10,
        plot(theta,w,'-k','Linewidth',3);
    end;
    
    if ijk==10,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,2*width,col);
    end;

    
    if ijk==5,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==15,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
    if ijk==20,
        ind=find(theta>0);
        ind=min(ind);
        x=theta(ind); y=w(ind);
        v=[theta(ind+1)-theta(ind); w(ind+1)-w(ind)];
        arrow(A,B,C,D,x,y,v,epsilon,width,col);
    end;
    
end;

hold off;

shg;
