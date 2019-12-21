clear; clf;
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

i_ext=1.2; 
    
N=10;               % number of neurons in the network

t_final=30;
dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

initial_vector=rtm_init(i_ext,[N-1:-1:0]/N+1/(2*N));

z=zeros(N,m_steps+1);
v=z;
m=z;
h=z;
n=z;

v(:,1)=initial_vector(:,1);
m(:,1)=m_inf(v(:,1));
h(:,1)=initial_vector(:,2);
n(:,1)=initial_vector(:,3);


for k=1:m_steps,
    
    v_inc=(g_k*n(:,k).^4.*(v_k-v(:,k))+g_na*m(:,k).^3.*h(:,k).*(v_na-v(:,k))+...
        g_l*(v_l-v(:,k))+i_ext)/c;
    n_inc=alpha_n(v(:,k)).*(1-n(:,k))-beta_n(v(:,k)).*n(:,k);
    h_inc=alpha_h(v(:,k)).*(1-h(:,k))-beta_h(v(:,k)).*h(:,k);
    
    v_tmp=v(:,k)+dt05*v_inc;
    h_tmp=h(:,k)+dt05*h_inc;
    n_tmp=n(:,k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
        g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v(:,k+1)=v(:,k)+dt*v_inc;
    h(:,k+1)=h(:,k)+dt*h_inc;
    n(:,k+1)=n(:,k)+dt*n_inc;
    m(:,k+1)=m_inf(v(:,k+1));
    
end;

t=[0:m_steps]*dt;
subplot(311);
set(gca,'Fontsize',14);

% Plot voltage traces:

for k=1:N,
    plot(t,v(k,:),'-k','Linewidth',1);
    hold on;
end;
hold off;
axis([0,t_final,-100,50]);
shg;
ylabel('$v$ [mV]','Fontsize',18);


g_syn=.30; v_syn=-75; tau_syn=9;

initial_vector=rtm_init(i_ext,[N-1:-1:0]/N+1/(2*N));

z=zeros(N,m_steps+1);
v=z;
m=z;
h=z;
n=z;

v(:,1)=initial_vector(:,1);
m(:,1)=m_inf(v(:,1));
h(:,1)=initial_vector(:,2);
n(:,1)=initial_vector(:,3);


for k=1:m_steps,
    
    v_inc=(g_k*n(:,k).^4.*(v_k-v(:,k))+g_na*m(:,k).^3.*h(:,k).*(v_na-v(:,k))+...
        g_l*(v_l-v(:,k))+i_ext+g_syn*exp(-(k-1)*dt/tau_syn)*(v_syn-v(:,k)))/c;
    n_inc=alpha_n(v(:,k)).*(1-n(:,k))-beta_n(v(:,k)).*n(:,k);
    h_inc=alpha_h(v(:,k)).*(1-h(:,k))-beta_h(v(:,k)).*h(:,k);
    
    v_tmp=v(:,k)+dt05*v_inc;
    h_tmp=h(:,k)+dt05*h_inc;
    n_tmp=n(:,k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
        g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext+g_syn*exp(-(k-1/2)*dt/tau_syn)*(v_syn-v_tmp))/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v(:,k+1)=v(:,k)+dt*v_inc;
    h(:,k+1)=h(:,k)+dt*h_inc;
    n(:,k+1)=n(:,k)+dt*n_inc;
    m(:,k+1)=m_inf(v(:,k+1));
    
end;

t=[0:m_steps]*dt;
subplot(312);
set(gca,'Fontsize',14);

% Plot spike rastergram:

for k=1:N,
    plot(t,v(k,:),'-k','Linewidth',1);
    hold on;
end;
hold off;
axis([0,t_final,-100,50]);
shg;
ylabel('$v$ [mV]','Fontsize',18);
    


g_syn=2.25; v_syn=-75; tau_syn=1;

initial_vector=rtm_init(i_ext,[N-1:-1:0]/N+1/(2*N));

z=zeros(N,m_steps+1);
v=z;
m=z;
h=z;
n=z;

v(:,1)=initial_vector(:,1);
m(:,1)=m_inf(v(:,1));
h(:,1)=initial_vector(:,2);
n(:,1)=initial_vector(:,3);


for k=1:m_steps,
    
    v_inc=(g_k*n(:,k).^4.*(v_k-v(:,k))+g_na*m(:,k).^3.*h(:,k).*(v_na-v(:,k))+...
        g_l*(v_l-v(:,k))+i_ext+g_syn*exp(-(k-1)*dt/tau_syn)*(v_syn-v(:,k)))/c;
    n_inc=alpha_n(v(:,k)).*(1-n(:,k))-beta_n(v(:,k)).*n(:,k);
    h_inc=alpha_h(v(:,k)).*(1-h(:,k))-beta_h(v(:,k)).*h(:,k);
    
    v_tmp=v(:,k)+dt05*v_inc;
    h_tmp=h(:,k)+dt05*h_inc;
    n_tmp=n(:,k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
        g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext+g_syn*exp(-(k-1/2)*dt/tau_syn)*(v_syn-v_tmp))/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v(:,k+1)=v(:,k)+dt*v_inc;
    h(:,k+1)=h(:,k)+dt*h_inc;
    n(:,k+1)=n(:,k)+dt*n_inc;
    m(:,k+1)=m_inf(v(:,k+1));
    
end;


t=[0:m_steps]*dt;
subplot(313);
set(gca,'Fontsize',14);

% Plot spike rastergram:

for k=1:N,
    plot(t,v(k,:),'-k','Linewidth',1);
    hold on;
end;
hold off;
axis([0,t_final,-100,50]);
shg;
xlabel('$t$ [ms]','Fontsize',18); 
ylabel('$v$ [mV]','Fontsize',18);
    
    
    
