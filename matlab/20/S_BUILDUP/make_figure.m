clear; clf;

tau_r=10; 
tau_d=300;
% tau_peak=20;
% tau_dq=tau_d_q_function(tau_d,tau_r,tau_peak);
tau_dq=5;

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;
i_ext=0.2; 


t_final=2000;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; s=z; q=z;

v(1)=-70;
m(1)=m_inf(v(1));
h(1)=h_inf(v(1));
n(1)=n_inf(v(1));
s(1)=0;
q(1)=0;


num_spikes=0;
for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    q_inc=(1+tanh(v(k)/10))/2*(1-q(k))/0.1-q(k)/tau_dq;
    s_inc=q(k)*(1-s(k))/tau_r-s(k)/tau_d;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    q_tmp=q(k)+dt05*q_inc;
    s_tmp=s(k)+dt05*s_inc;
    
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    q_inc=(1+tanh(v_tmp/10))/2*(1-q_tmp)/0.1-q_tmp/tau_dq;
    s_inc=q_tmp*(1-s_tmp)/tau_r-s_tmp/tau_d;
    
    
    v(k+1)=v(k)+dt*v_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    m(k+1)=m_inf(v(k+1));
    q(k+1)=q(k)+dt*q_inc;
    s(k+1)=s(k)+dt*s_inc;
    
    if v(k+1)<=-20 & v(k)>-20,
        num_spikes=num_spikes+1;
        ts=(k-1)*dt*(-20-v(k+1))+k*dt*(v(k)+20);
        ts=ts/(v(k)-v(k+1));
        t_spikes(num_spikes)=ts;
    end;
    
end;

period=t_spikes(num_spikes)-t_spikes(num_spikes-1)


t=[0:m_steps]*dt;
subplot(211);
set(gca,'Fontsize',16);
plot(t,v,'-k','Linewidth',2);
ylabel('$v$ [mV]','Fontsize',20);
axis([0,t_final,-100,100]);
    
subplot(212);
set(gca,'Fontsize',16);
plot(t,s,'-k','Linewidth',2);
ylabel('$s$','Fontsize',20);
axis([0,t_final,0,1]);
xlabel('$t$ [ms]','Fontsize',20);
shg;

