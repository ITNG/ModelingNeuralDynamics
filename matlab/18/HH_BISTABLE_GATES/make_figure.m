clear;
c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;

i_ext=8;

t_final=40;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);

v_vec=-100+[0:100000]/100000*150;
f=@(v) g_na*m_inf(v).^3.*h_inf(v).*(v_na-v)+g_k*n_inf(v).^4.*(v_k-v)+ ...
    g_l*(v_l-v)+i_ext;

v_left=-100;
v_right=50;
while v_right-v_left>10^(-12),
    v_c=(v_left+v_right)/2;
    if f(v_c)*f(v_right)<0,
        v_left=v_c;
    else
        v_right=v_c;
    end;
end;
v_star=(v_left+v_right)/2;
m_star=m_inf(v_star);
h_star=h_inf(v_star);
n_star=n_inf(v_star);

v(1)=v_star+5;
m(1)=m_star;
h(1)=h_star; 
n(1)=n_star;

for k=1:m_steps,
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    m_inc=alpha_m(v(k))*(1-m(k))-beta_m(v(k))*m(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m(k)+dt05*m_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m(k)+dt*m_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    
end;

t=[0:m_steps]*dt;

subplot(311);
plot(t,m,'-k','Linewidth',2);
hold on;
plot(t,m_inf(v),'--b','Linewidth',2);
plot([0,t_final],[m_star,m_star],'-r','Linewidth',2);
hold off;
set(gca,'Fontsize',14);
ylabel('$m$','Fontsize',18);
axis([0,t_final,0,1]);
shg;

subplot(312);
plot(t,h,'-k','Linewidth',2);
hold on;
plot(t,h_inf(v),'--b','Linewidth',2);
plot([0,t_final],[h_star,h_star],'-r','Linewidth',2);
ind=find(h>h_star);
for i=1:length(ind);
    plot(t(ind(i)),0,'.m','Markersize',20);
end;
hold off;
set(gca,'Fontsize',14);
ylabel('$h$','Fontsize',18);
axis([0,t_final,0,1]);

subplot(313);
plot(t,n,'-k','Linewidth',2);
hold on;
plot(t,n_inf(v),'--b','Linewidth',2);
plot([0,t_final],[n_star,n_star],'-r','Linewidth',2);
ind=find(n<n_star);
for i=1:length(ind);
    plot(t(ind(i)),0,'.m','Markersize',20);
end;
hold off;
set(gca,'Fontsize',14);
ylabel('$n$','Fontsize',18);
xlabel('$t$ [ms]','Fontsize',18);
axis([0,t_final,0,1]);

shg;

shg;
    
    
