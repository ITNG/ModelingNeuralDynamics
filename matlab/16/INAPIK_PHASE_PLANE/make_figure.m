clear; clf;

c=1;
g_na=20;
g_k=10; 
g_l=8;
v_na=60;
v_k=-90;
v_l=-80;
tau_n=0.15;


i_ext=-5;

t_final=2.5;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
z=zeros(m_steps+1,1);
v=z; m=z; n=z;


v(1)=-50; 
m(1)=m_inf(v(1));
n(1)=0;

for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
end;

subplot(221);
plot(v,n,'-k','Linewidth',2);
shg; hold on;

f=@(v) g_na*m_inf(v).*(v_na-v)+g_k*n_inf(v).*(v_k-v)+g_l*(v_l-v);

w=-100+[0:10000]/10000*150;

k=[1:length(w)-1];
ind=find((f(w(k))+i_ext).*(f(w(k+1))+i_ext)<=0);
for klm=1:length(ind),
    j=ind(klm);
    w_low=w(j); w_high=w(j+1);
    while w_high-w_low>10^(-12),
        w_c=(w_low+w_high)/2;
        if (f(w_c)+i_ext)*(f(w_high)+i_ext)<=0,
            w_low=w_c;
        else
            w_high=w_c;
        end;
    end;
    v_c=(w_low+w_high)/2;
    n_c=n_inf(v_c);
    J=zeros(2,2);
    J(1,1)=g_na*m_inf_p(v_c)*(v_na-v_c)-g_na*m_inf(v_c)-g_k*n_c-g_l;
    J(1,2)=g_k*(v_k-v_c);
    J(2,1)=n_inf_p(v_c)/tau_n;
    J(2,2)=-1/tau_n;
    E=eig(J);
    if abs(imag(E(1)))<10^(-12) & real(E(1))<0 & real(E(2))<0,
        % stable node
        plot(v_c,n_c,'.k','Markersize',30);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))>0 & real(E(2))>0,
        % unstable node
        plot(v_c,n_c,'.b','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))<0,
        % stable spiral
        plot(v_c,n_c,'.r','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))>0,
        % unstable spiral
        plot(v_c,n_c,'og','Markersize',8);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))*real(E(2))<0,
        % saddle
        plot(v_c,n_c,'om','Markersize',8);
    end;

end;

A=-75; B=0; C=-0.15; D=0.85;
axis([A,B,C,D]);

t0=0.4;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
epsilon=0.07;
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

t0=0.7;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

t0=1.8;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
theta=0.002;
v0=([cos(theta) -sin(theta);
    sin(theta) cos(theta)]*v0')';
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

hold off;

axis('square');
set(gca,'Fontsize',16);
%xlabel('$v$','Fontsize',20);
set(gca,'Xtick',[]);
ylabel('$n$','Fontsize',20);
title('$I =-5$','Fontsize',20);

t_final=3;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
z=zeros(m_steps+1,1);
v=z; m=z; n=z;

i_ext=-1.4;

v(1)=-52.5; 
m(1)=m_inf(v(1));
n(1)=0; 

for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
end;

subplot(222);                          
v_plot=v;                             
n_plot=n;
plot(v_plot,n_plot,'-k','Linewidth',2);
shg; hold on;

f=@(v) g_na*m_inf(v).*(v_na-v)+g_k*n_inf(v).*(v_k-v)+g_l*(v_l-v);

w=-100+[0:10000]/10000*150;

k=[1:length(w)-1];
ind=find((f(w(k))+i_ext).*(f(w(k+1))+i_ext)<=0);
for klm=1:length(ind),
    j=ind(klm);
    w_low=w(j); w_high=w(j+1);
    while w_high-w_low>10^(-12),
        w_c=(w_low+w_high)/2;
        if (f(w_c)+i_ext)*(f(w_high)+i_ext)<=0,
            w_low=w_c;
        else
            w_high=w_c;
        end;
    end;
    v_c=(w_low+w_high)/2;
    n_c=n_inf(v_c);
    J=zeros(2,2);
    J(1,1)=g_na*m_inf_p(v_c)*(v_na-v_c)-g_na*m_inf(v_c)-g_k*n_c-g_l;
    J(1,2)=g_k*(v_k-v_c);
    J(2,1)=n_inf_p(v_c)/tau_n;
    J(2,2)=-1/tau_n;
    E=eig(J);
    if abs(imag(E(1)))<10^(-12) & real(E(1))<0 & real(E(2))<0,
        % stable node
        plot(v_c,n_c,'.k','Markersize',30);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))>0 & real(E(2))>0,
        % unstable node
        plot(v_c,n_c,'.b','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))<0,
        % stable spiral
        plot(v_c,n_c,'.r','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))>0,
        % unstable spiral
        plot(v_c,n_c,'og','Markersize',8);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))*real(E(2))<0,
        % saddle
        plot(v_c,n_c,'om','Markersize',8);
    end;

end;

A=-75; B=0; C=-0.15; D=0.85;
axis([A,B,C,D]);

t0=0.5;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

t0=0.9;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

hold off;

axis('square');
set(gca,'Fontsize',16);
%xlabel('$v$','Fontsize',20);
%ylabel('$n$','Fontsize',20);
set(gca,'Ytick',[]);
set(gca,'Xtick',[]);
title('$I =-1.4$','Fontsize',20);



t_final=20;
dt=0.001; dt05=dt/2;
m_steps=round(t_final/dt);
z=zeros(m_steps+1,1);
v=z; m=z; n=z;

i_ext=2;

v(1)=-30; 
m(1)=m_inf(v(1));
n(1)=0.2; 

for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
end;

subplot(223);
ind=[round((t_final-1.2)/dt):m_steps+1]; % We display only the last
v_plot=v(ind);                           % 1.2 ms of the trajectory. By
n_plot=n(ind);                           % trial and error, we found that
plot(v_plot,n_plot,'-k','Linewidth',2);  % that's enough for one complete
shg; hold on;                            % pass through the limit cycle,                                        
                                         % but not much more.
                                         
f=@(v) g_na*m_inf(v).*(v_na-v)+g_k*n_inf(v).*(v_k-v)+g_l*(v_l-v);

w=-100+[0:10000]/10000*150;

k=[1:length(w)-1];
ind=find((f(w(k))+i_ext).*(f(w(k+1))+i_ext)<=0);
for klm=1:length(ind),
    j=ind(klm);
    w_low=w(j); w_high=w(j+1);
    while w_high-w_low>10^(-12),
        w_c=(w_low+w_high)/2;
        if (f(w_c)+i_ext)*(f(w_high)+i_ext)<=0,
            w_low=w_c;
        else
            w_high=w_c;
        end;
    end;
    v_c=(w_low+w_high)/2;
    n_c=n_inf(v_c);
    J=zeros(2,2);
    J(1,1)=g_na*m_inf_p(v_c)*(v_na-v_c)-g_na*m_inf(v_c)-g_k*n_c-g_l;
    J(1,2)=g_k*(v_k-v_c);
    J(2,1)=n_inf_p(v_c)/tau_n;
    J(2,2)=-1/tau_n;
    E=eig(J);
    if abs(imag(E(1)))<10^(-12) & real(E(1))<0 & real(E(2))<0,
        % stable node
        plot(v_c,n_c,'.k','Markersize',30);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))>0 & real(E(2))>0,
        % unstable node
        plot(v_c,n_c,'.b','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))<0,
        % stable spiral
        plot(v_c,n_c,'.r','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))>0,
        % unstable spiral
        plot(v_c,n_c,'og','Markersize',8);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))*real(E(2))<0,
        % saddle
        plot(v_c,n_c,'om','Markersize',8);
    end;

end;

% Here we put two arrows into the computed limit cycle. 
% By trial and error, we found ways of placing the arrows
% that make nice plots:

A=-75; B=0; C=-0.15; D=0.85;
axis([A,B,C,D]);
t0=t_final-1.6;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

t0=t_final-1.9;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

hold off;

axis('square');
set(gca,'Fontsize',16);
xlabel('$v$','Fontsize',20);
ylabel('$n$','Fontsize',20);
title('$I=2$','Fontsize',20);





i_ext=4.4;

v(1)=-30; 
m(1)=m_inf(v(1));
n(1)=0.2; 

for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
end;

subplot(224);
ind=[round((t_final-2)/dt):m_steps+1];
v_plot=v(ind);
n_plot=n(ind);
plot(v_plot,n_plot,'-k','Linewidth',2);
shg;
hold on;

f=@(v) g_na*m_inf(v).*(v_na-v)+g_k*n_inf(v).*(v_k-v)+g_l*(v_l-v);

w=-100+[0:10000]/10000*150;

k=[1:length(w)-1];
ind=find((f(w(k))+i_ext).*(f(w(k+1))+i_ext)<=0);
for klm=1:length(ind),
    j=ind(klm);
    w_low=w(j); w_high=w(j+1);
    while w_high-w_low>10^(-12),
        w_c=(w_low+w_high)/2;
        if (f(w_c)+i_ext)*(f(w_high)+i_ext)<=0,
            w_low=w_c;
        else
            w_high=w_c;
        end;
    end;
    v_c=(w_low+w_high)/2;
    n_c=n_inf(v_c);
    J=zeros(2,2);
    J(1,1)=g_na*m_inf_p(v_c)*(v_na-v_c)-g_na*m_inf(v_c)-g_k*n_c-g_l;
    J(1,2)=g_k*(v_k-v_c);
    J(2,1)=n_inf_p(v_c)/tau_n;
    J(2,2)=-1/tau_n;
    E=eig(J);
    if abs(imag(E(1)))<10^(-12) & real(E(1))<0 & real(E(2))<0,
        plot(v_c,n_c,'.k','Markersize',30);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))>0 & real(E(2))>0,
        plot(v_c,n_c,'.b','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))<0,
        plot(v_c,n_c,'.r','Markersize',30);
    end;
    if abs(imag(E(1)))>10^(-12) & real(E(1))>0,
        plot(v_c,n_c,'og','Markersize',8);
    end;
    if abs(imag(E(1)))<10^(-12) & real(E(1))*real(E(2))<0,
        plot(v_c,n_c,'om','Markersize',8);
    end;

end;
A=-75; B=0; C=-0.15; D=0.85;
axis([A,B,C,D]);

t0=t_final-1;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);

t0=t_final-1.3;
k0=round(t0/dt)+1;
x0=v(k0);
y0=n(k0);
v0(1)=(v(k0+1)-v(k0-1))/(2*dt);
v0(2)=(n(k0+1)-n(k0-1))/(2*dt);
width=2;
col='k';
arrow(A,B,C,D,x0,y0,v0,epsilon,width,col);


hold off;

axis('square');
set(gca,'Fontsize',16);
xlabel('$v$','Fontsize',20);
set(gca,'Ytick',[]);
title('$I=4.4$','Fontsize',20);



