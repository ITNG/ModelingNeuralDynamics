

clear; clf;
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;



t_final=30;
dt=0.005;
dt05=dt/2;
m_steps=round(t_final/dt);

i_ext=1; 

v(1)=-75.582445796204553;
m(1)=m_inf(v(1));
n(1)=0.005935897840905;
h(1)=1-n(1);
A=-105; B=55; C=-0.05; D=0.8;

subplot(231);
hold on;
for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    h_tmp=1-n_tmp;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
    h(k+1)=1-n(k+1);
    m(k+1)=m_inf(v(k+1));
    if v(k)<0&v(k+1)>0,
        x=v(k);
        y=n(k);
        vv=[v(k+1)-v(k); n(k+1)-n(k)];
        epsilon=0.1;
        width=2;
        arrow(A,B,C,D,x,y,vv,epsilon,width,'-k')
    end;
    if v(k)>-30&v(k+1)<-30,
        x=v(k);
        y=n(k);
        vv=[v(k+1)-v(k); n(k+1)-n(k)];
        epsilon=0.1;
        width=2;
        arrow(A,B,C,D,x,y,vv,epsilon,width,'-k')
    end;
end;

plot(v,n,'-k','Linewidth',2);


shg; 
hold off;
set(gca,'Fontsize',12);
set(gca,'Box','on');
xlabel('$v$ [mV]','Fontsize',16); ylabel('$n$','Fontsize',16);
axis([A,B,C,D]);
axis('square');



i_ext=0;


f=@(v) g_k*n_inf(v).^4.*(v_k-v)+g_na*m_inf(v).^3.*h_inf(v).*(v_na-v) ...
    +g_l*(v_l-v)+i_ext;
v_left=-63; v_right=-62;
while v_right-v_left>10^(-14),
    v_c=(v_left+v_right)/2;
    if f(v_c)*f(v_left)>0,
        v_left=v_c;
    else
        v_right=v_c;
    end;
end;
v_star=v_c;
n_star=n_inf(v_c);

v_left=-75; v_right=-65;
while v_right-v_left>10^(-12),
    v_c=(v_left+v_right)/2;
    if f(v_c)*f(v_left)>0,
        v_left=v_c;
    else
        v_right=v_c;
    end;
end;
v_0=v_c;
n_0=n_inf(v_c);


t_final=300;
dt=0.005;
dt05=dt/2;
m_steps=round(t_final/dt);

v(1)=v_star+.5;
m(1)=m_inf(v(1));
n(1)=n_star+0.005;
h(1)=1-n(1);

subplot(232);
hold on;

for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    h_tmp=1-n_tmp;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
    h(k+1)=1-n(k+1);
    m(k+1)=m_inf(v(k+1));
    
    if v(k)<0&v(k+1)>0,
        x=v(k);
        y=n(k);
        vv=[v(k+1)-v(k); n(k+1)-n(k)];
        epsilon=0.1;
        width=2;
        arrow(A,B,C,D,x,y,vv,epsilon,width,'-k')
    end;
    if v(k)>-30&v(k+1)<-30,
        x=v(k);
        y=n(k);
        vv=[v(k+1)-v(k); n(k+1)-n(k)];
        epsilon=0.1;
        width=2;
        arrow(A,B,C,D,x,y,vv,epsilon,width,'-k')
    end;
    
end;



plot(v,n,'-k','Linewidth',2);
hold off;
set(gca,'Fontsize',12);
set(gca,'Box','on');
xlabel('$v$ [mV]','Fontsize',16); %ylabel('$n$','Fontsize',16);
axis([A,B,C,D]);
axis('square');

subplot(233);
A=-75; B=-50; C=0; D=0.15;

plot(v,n,'-k','Linewidth',2);
hold on;
kc=[1:m_steps];
kr=[2:m_steps+1];
ind=find(v(kc)<-55&v(kr)>=-55);
k=ind(1);
x=v(k);
y=n(k);
vv=[v(k+1)-v(k); n(k+1)-n(k)];
epsilon=0.1;
width=2;
arrow(A,B,C,D,x,y,vv,epsilon,width,'-k')
ind=find(v(kc)<-70&v(kr)>=-70);
k=ind(1);
x=v(k);
y=n(k);
vv=[v(k+1)-v(k); n(k+1)-n(k)];
epsilon=0.1;
width=2;
arrow(A,B,C,D,x,y,vv,epsilon,width,'-k')





v(1)=v_star-.5;
m(1)=m_inf(v(1));
n(1)=n_star-0.005;
h(1)=1-n(1);


for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    h_tmp=1-n_tmp;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
    h(k+1)=1-n(k+1);
    m(k+1)=m_inf(v(k+1));
    
    if v(k)>-65&v(k+1)<-65,
        x=v(k);
        y=n(k);
        vv=[v(k+1)-v(k); n(k+1)-n(k)];
        epsilon=0.1;
        width=2;
        arrow(A,B,C,D,x,y,vv,epsilon,width,'-r')
    end;
    
end;



hold on;
plot(v,n,'-r','Linewidth',2);
plot(v_star,n_star,'ok','Markersize',8,'Linewidth',1);
plot(v_0,n_0,'ok','Markersize',8,'MarkerFaceColor','k');
shg;
hold off;
set(gca,'Fontsize',12);
xlabel('$v$ [mV]','Fontsize',16); %ylabel('$n$','Fontsize',16);
axis([A,B,C,D]);
axis('square');





