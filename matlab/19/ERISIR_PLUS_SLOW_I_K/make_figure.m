clear; clf;

c=1;
g_k=224;
g_na=112;
g_l=0.5;
v_k=-90;
v_na=60;
v_l=-70;

i_ext=7.5;

g_k_slow=1.5; tau_n_slow=100;


t_final=1000;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; n_slow=z;

v(1)=-70;
m(1)=m_inf(v(1));
h(1)=h_inf(v(1));
n(1)=n_inf(v(1));
n_slow(1)=n_slow_inf(v(1));

for k=1:m_steps,
    
    v_inc=(g_k*n(k)^2*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
           g_k_slow*n_slow(k)*(v_k-v(k))+ ...
           g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_slow_inc=(n_slow_inf(v(k))-n_slow(k))/tau_n_slow;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    n_slow_tmp=n_slow(k)+dt05*n_slow_inc;
    
    v_inc=(g_k*n_tmp^2*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_k_slow*n_slow_tmp*(v_k-v_tmp)+ ...
           g_l*(v_l-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    n_slow_inc=(n_slow_inf(v_tmp)-n_slow_tmp)/tau_n_slow;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    n_slow(k+1)=n_slow(k)+dt*n_slow_inc;
    
end;

t=[0:m_steps]*dt;
subplot(211);
plot(t,v,'-k','Linewidth',2);
axis([0,t_final,-95,55]);
shg;
hold off;
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$v$ [mV]','Fontsize',20);
    
