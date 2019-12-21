clear; clf;
c=1;                    % parameters for RTM neuron
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

i_ext=0.3; 

t_final=300;
dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

tau_r=0.5; tau_peak=0.5; tau_d=2;       % AMPA-like synaptic input pulse
tau_d_q=tau_d_q_function(tau_d,tau_r,tau_peak);
g_syn=0.1;
initial_vector=rtm_init(i_ext,0.2);
                                
z=zeros(m_steps+1,1); v=z; m=z; h=z; n=z; q=z; s=z;

v(1)=initial_vector(1,1);
m(1)=m_inf(v(1));
h(1)=initial_vector(1,2);
n(1)=initial_vector(1,3);

for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+ ...
           g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
           g_l*(v_l-v(k))+ ...
           i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+ ...
           g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+ ...
           i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;

end;

t=[0:m_steps]*dt;
subplot(211);
set(gca,'Fontsize',16);
plot(t,v,'-b','Linewidth',4);
ylabel('$v$ [mV]','Fontsize',20);
    
    
initial_vector=rtm_init(i_ext,0.2);
                                
v(1)=initial_vector(1,1);
m(1)=m_inf(v(1));
h(1)=initial_vector(1,2);
n(1)=initial_vector(1,3);
q(1)=0;                     % The synaptic input comes in later. 
s(1)=0;


for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+ ...
           g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
           g_l*(v_l-v(k))+ ...
           g_syn*s(k)*(-v(k)) + ...
           i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    q_inc=-q(k)/tau_d_q;
    s_inc=q(k)*(1-s(k))/tau_r-s(k)/tau_d;
    
    v_tmp=v(k)+dt05*v_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    q_tmp=q(k)+dt05*q_inc;
    s_tmp=s(k)+dt05*s_inc;
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+ ...
           g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+ ...
           g_syn*s_tmp*(-v_tmp)+ ...
           i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    q_inc=-q_tmp/tau_d_q;
    s_inc=q_tmp*(1-s_tmp)/tau_r-s_tmp/tau_d;
    
    v(k+1)=v(k)+dt*v_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    m(k+1)=m_inf(v(k+1));
    q(k+1)=q(k)+dt*q_inc;
    s(k+1)=s(k)+dt*s_inc;
    
    if abs((k+1)*dt-100)<10^(-6)
        q(k+1)=1;       % Here, at time t=100, comes the synaptic
                        % input pulse. 
    end;
    
end;

hold on;
plot([100,100],[-100,50],'--k','Linewidth',1);
plot(t,v,'-r','Linewidth',1);
hold off;


initial_vector=rtm_init(i_ext,0.2);
                                
z=zeros(m_steps+1,1); v=z; m=z; h=z; n=z; q=z; s=z;

v(1)=initial_vector(1,1);
m(1)=m_inf(v(1));
h(1)=initial_vector(1,2);
n(1)=initial_vector(1,3);



for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+ ...
           g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
           g_l*(v_l-v(k))+ ...
           i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+ ...
           g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+ ...
           i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;

end;

t=[0:m_steps]*dt;
subplot(212);
set(gca,'Fontsize',16);
plot(t,v,'-b','Linewidth',4);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$v$ [mV]','Fontsize',20);
    
    
initial_vector=rtm_init(i_ext,0.2);
                                
v(1)=initial_vector(1,1);
m(1)=m_inf(v(1));
h(1)=initial_vector(1,2);
n(1)=initial_vector(1,3);
q(1)=0;                     % The synaptic input pulse comes in later.
s(1)=0;


for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+ ...
           g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
           g_l*(v_l-v(k))+ ...
           g_syn*s(k)*(-v(k)) + ...
           i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    q_inc=-q(k)/tau_d_q;
    s_inc=q(k)*(1-s(k))/tau_r-s(k)/tau_d;
    
    v_tmp=v(k)+dt05*v_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    q_tmp=q(k)+dt05*q_inc;
    s_tmp=s(k)+dt05*s_inc;
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+ ...
           g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+ ...
           g_syn*s_tmp*(-v_tmp)+ ...
           i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    q_inc=-q_tmp/tau_d_q;
    s_inc=q_tmp*(1-s_tmp)/tau_r-s_tmp/tau_d;
    
    v(k+1)=v(k)+dt*v_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    m(k+1)=m_inf(v(k+1));
    q(k+1)=q(k)+dt*q_inc;
    s(k+1)=s(k)+dt*s_inc;
    
    if abs((k+1)*dt-120)<10^(-6)    % This time the synaptic input pulse
        q(k+1)=1;                   % arrives at time t=120. 
    end;
    
end;

hold on;
plot([120,120],[-100,50],'--k','Linewidth',1);
plot(t,v,'-r','Linewidth',1);
hold off;
shg; 
    
