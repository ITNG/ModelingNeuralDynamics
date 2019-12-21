clear; clf;

c=1;
g_na=20;
g_k=10; 
g_l=8;
v_na=60;
v_k=-90;
v_l=-80;
tau_n=0.15;

i_ext=7;

g_k_slow=4; tau_n_slow=20;

t_final=100; dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); v=z; m=z; n=z; n_slow=z;
v(1)=-70; m(1)=m_inf(v(1)); n(1)=0.6; n_slow(1)=0;

for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+ ...
        g_k_slow*n_slow(k)*(v_k-v(k))+ ...
        g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    n_slow_inc=(n_slow_inf(v(k))-n_slow(k))/tau_n_slow;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    n_slow_tmp=n_slow(k)+dt05*n_slow_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+ ...
        g_k_slow*n_slow_tmp*(v_k-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    n_slow_inc=(n_slow_inf(v_tmp)-n_slow_tmp)/tau_n_slow;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
    n_slow(k+1)=n_slow(k)+dt*n_slow_inc;
    
end;

subplot(211); t=[0:m_steps]*dt;
plot(t,v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$v$ [mV]','Fontsize',20);
shg;