% The trajectories plotted here, and the locations of the arrows, were 
% put in with a lot of trial and error to make the figure look nice.
% Don't try to detect a system, there was none. 


tau_m=1; I=0.3; tau_I=9;


subplot(111);

M=6;
for ijk=1:M, 
    k=1;
    clear theta;
    clear g_syn;
    theta(1)=-pi;
    g_syn(1)=ijk/M;
    dt=0.01;
    dt05=dt/2;

    while theta(k)>=-pi&theta(k)<=pi,
        theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
            -g_syn(k)*sin(theta(k));
        g_syn_inc=-g_syn(k)/tau_I;
        theta_tmp=theta(k)+dt05*theta_inc;
        g_syn_tmp=g_syn(k)+dt05*g_syn_inc;
        theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
            -g_syn_tmp*sin(theta_tmp);
        g_syn_inc=-g_syn_tmp/tau_I;
        theta(k+1)=theta(k)+dt*theta_inc;
        g_syn(k+1)=g_syn(k)+dt*g_syn_inc;
        k=k+1;
    end;

    N=k;

    plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);
    hold on;
end;

M=6;
for ijk=1:M,
    k=1;
    clear theta;
    clear g_syn;
    theta(1)=pi;
    g_syn(1)=ijk/M;
    dt=0.01;
    dt05=dt/2;

    while theta(k)>=-pi&theta(k)<=pi,
        theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
            -g_syn(k)*sin(theta(k));
        g_syn_inc=-g_syn(k)/tau_I;
        theta_tmp=theta(k)-dt05*theta_inc;  % Notice the minus time. We plot
                                            % trajectories going backwards
                                            % in time in this loop. 
        g_syn_tmp=g_syn(k)-dt05*g_syn_inc;
        theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
            -g_syn_tmp*sin(theta_tmp);
        g_syn_inc=-g_syn_tmp/tau_I;
        theta(k+1)=theta(k)-dt*theta_inc;
        g_syn(k+1)=g_syn(k)-dt*g_syn_inc;
        k=k+1;
    end;

    N=k;

    plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);

end;

M=6;
for ijk=1:M,
    k=1;
    clear theta;
    clear g_syn;
    theta(1)=0;
    g_syn(1)=ijk/(M+1);
    dt=0.01;
    dt05=dt/2;

    while theta(k)>=-pi&theta(k)<=pi,
        theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
            -g_syn(k)*sin(theta(k));
        g_syn_inc=-g_syn(k)/tau_I;
        theta_tmp=theta(k)+dt05*theta_inc;
        g_syn_tmp=g_syn(k)+dt05*g_syn_inc;
        theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
            -g_syn_tmp*sin(theta_tmp);
        g_syn_inc=-g_syn_tmp/tau_I;
        theta(k+1)=theta(k)+dt*theta_inc;
        g_syn(k+1)=g_syn(k)+dt*g_syn_inc;
        k=k+1;
    end;

    N=k;

    plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);
end;
for ijk=1:M,
    k=1;
    clear theta;
    clear g_syn;
    theta(1)=0;
    g_syn(1)=ijk/(M+1);
    dt=0.01;
    dt05=dt/2;

    while theta(k)>=-pi&theta(k)<=pi,
        theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
            -g_syn(k)*sin(theta(k));
        g_syn_inc=-g_syn(k)/tau_I;
        theta_tmp=theta(k)-dt05*theta_inc;  % About the minus sign, see above.
        g_syn_tmp=g_syn(k)-dt05*g_syn_inc;
        theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
            -g_syn_tmp*sin(theta_tmp);
        g_syn_inc=-g_syn_tmp/tau_I;
        theta(k+1)=theta(k)-dt*theta_inc;
        g_syn(k+1)=g_syn(k)-dt*g_syn_inc;
        k=k+1;
    end;

    N=k;

    plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);
end;


k=1;
clear theta;
clear g_syn;
theta(1)=-2;
g_syn(1)=1;
dt=0.01;
dt05=dt/2;

while theta(k)>=-pi&theta(k)<=pi,
    theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
        -g_syn(k)*sin(theta(k));
    g_syn_inc=-g_syn(k)/tau_I;
    theta_tmp=theta(k)+dt05*theta_inc;
    g_syn_tmp=g_syn(k)+dt05*g_syn_inc;
    theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
        -g_syn_tmp*sin(theta_tmp);
    g_syn_inc=-g_syn_tmp/tau_I;
    theta(k+1)=theta(k)+dt*theta_inc;
    g_syn(k+1)=g_syn(k)+dt*g_syn_inc;
    k=k+1;
end;
N=k;
plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);
hold on;



k=1;
clear theta;
clear g_syn;
theta(1)=0;
g_syn(1)=1;
dt=0.01;
dt05=dt/2;

while theta(k)>=-pi&theta(k)<=pi,
    theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
        -g_syn(k)*sin(theta(k));
    g_syn_inc=-g_syn(k)/tau_I;
    theta_tmp=theta(k)+dt05*theta_inc;
    g_syn_tmp=g_syn(k)+dt05*g_syn_inc;
    theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
        -g_syn_tmp*sin(theta_tmp);
    g_syn_inc=-g_syn_tmp/tau_I;
    theta(k+1)=theta(k)+dt*theta_inc;
    g_syn(k+1)=g_syn(k)+dt*g_syn_inc;
    k=k+1;
end;
N=k;
plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);
hold on;


k=1;
clear theta;
clear g_syn;
theta(1)=pi;
g_syn(1)=0.075;
dt=0.01;
dt05=dt/2;

while theta(k)>=-pi&theta(k)<=pi,
    theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
        -g_syn(k)*sin(theta(k));
    g_syn_inc=-g_syn(k)/tau_I;
    theta_tmp=theta(k)-dt05*theta_inc;  % About the minus sign, see above.
    g_syn_tmp=g_syn(k)-dt05*g_syn_inc;
    theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
        -g_syn_tmp*sin(theta_tmp);
    g_syn_inc=-g_syn_tmp/tau_I;
    theta(k+1)=theta(k)-dt*theta_inc;
    g_syn(k+1)=g_syn(k)-dt*g_syn_inc;
    k=k+1;
end;
N=k;
plot(theta(1:N),g_syn(1:N),'-k','Linewidth',2);
hold on;


theta=0;
g_syn=4/7;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,2,'-k')


theta=-2;
g_syn=0.592;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,2,'-k')

theta=-2;
g_syn=0.145;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,2,'-k')

theta=0.5;
g_syn=0.203;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,2,'-k')

theta=2.5;
g_syn=0.18;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,2,'-k')

theta=2.85;
g_syn=0.52;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,2,'-k')


k=1;
clear theta;
clear g_syn;
theta(1)=-1.23;
g_syn(1)=1;
dt=0.01;
dt05=dt/2;

while theta(k)>=-pi&theta(k)<=pi,
    theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
        -g_syn(k)*sin(theta(k));
    g_syn_inc=-g_syn(k)/tau_I;
    theta_tmp=theta(k)+dt05*theta_inc;
    g_syn_tmp=g_syn(k)+dt05*g_syn_inc;
    theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
        -g_syn_tmp*sin(theta_tmp);
    g_syn_inc=-g_syn_tmp/tau_I;
    theta(k+1)=theta(k)+dt*theta_inc;
    g_syn(k+1)=g_syn(k)+dt*g_syn_inc;
    k=k+1;
end;
N=k;
plot(theta(1:N),g_syn(1:N),'-b','Linewidth',4);
g_ast=g_syn(N);
hold on;

k=1;
clear theta;
clear g_syn;
theta=-0.75;
g_syn=0.25;
theta_inc=-cos(theta)/tau_m+(2*I-g_syn)*(1+cos(theta)) ...
    -g_syn*sin(theta);
g_syn_inc=-g_syn/tau_I;
v=[theta_inc;g_syn_inc];
x=theta;
y=g_syn;
arrow(-pi,pi,0,1,x,y,v,0.05,4,'-b')

k=1;
clear theta;
clear g_syn;
theta(1)=1.25;
g_syn(1)=0.25;
dt=0.01;
dt05=dt/2;

while theta(k)>=-pi&theta(k)<=pi,
    theta_inc=-cos(theta(k))/tau_m+(2*I-g_syn(k))*(1+cos(theta(k))) ...
        -g_syn(k)*sin(theta(k));
    g_syn_inc=-g_syn(k)/tau_I;
    theta_tmp=theta(k)-dt05*theta_inc;  % About the minus sign, see above.
    g_syn_tmp=g_syn(k)-dt05*g_syn_inc;
    theta_inc=-cos(theta_tmp)/tau_m+(2*I-g_syn_tmp)*(1+cos(theta_tmp)) ...
        -g_syn_tmp*sin(theta_tmp);
    g_syn_inc=-g_syn_tmp/tau_I;
    theta(k+1)=theta(k)-dt*theta_inc;
    g_syn(k+1)=g_syn(k)-dt*g_syn_inc;
    k=k+1;
end;
N=k;
plot(theta(1:N),g_syn(1:N),'-r','Linewidth',4);

plot(pi,g_ast,'.k','Markersize',30);
text(pi+0.1,g_ast,'$g_\ast$','Fontsize',32);
hold off;

set(gca,'Fontsize',24);
axis([-pi,pi,0,1]); 
xlabel('$\theta$','Fontsize',32);
ylabel('$g_{\rm syn}$','Fontsize',32);
axis('square');
shg;
    

