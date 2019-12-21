clear; clf; rng('default'); rng(63806);

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

N=30;

% Set random external drives:
i_ext=0.25+rand(N,1)*0.1; 

tau_r=0.5; tau_peak=0.5; tau_d=2;
tau_d_q=tau_d_q_function(tau_d,tau_r,tau_peak);

% Set random synaptic onnection strengths:
g_syn=0.00625+rand(N,N)*0.0025;

% Remove autapses:
for i=1:N,
    g_syn(i,i)=0;
end;

t_final=200;
dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(N,1); v=z; h=z; n=z;
for i=1:N,
    initial_vector=rtm_init(i_ext(i),rand(1,1));
    v(i)=initial_vector(1);
    h(i)=initial_vector(2);
    n(i)=initial_vector(3);
end;

m=m_inf(v);
q=zeros(N,1);
s=zeros(N,1);

num_spikes=0;   % number of spikes so far
t_spikes=[];    % times of those spikes
i_spikes=[];    % neuronal indices associated with those
                % spikes (spike i was a spike of neuron
                % i_spikes(i))

for k=1:m_steps,
    
    v_inc=(g_k*n.^4.*(v_k-v)+g_na*m.^3.*h.*(v_na-v)+...
        g_l*(v_l-v)-(g_syn'*s).*v+i_ext)/c;
    h_inc=alpha_h(v).*(1-h)-beta_h(v).*h;
    n_inc=alpha_n(v).*(1-n)-beta_n(v).*n;
    q_inc=5*(1+tanh(v/10)).*(1-q)-q/tau_d_q;
    s_inc=q.*(1-s)/tau_r-s/tau_d;
    
    v_tmp=v+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h+dt05*h_inc;
    n_tmp=n+dt05*n_inc;
    q_tmp=q+dt05*q_inc;
    s_tmp=s+dt05*s_inc;
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
        g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)- ...
        (g_syn'*s_tmp).*v_tmp+i_ext)/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    q_inc=5*(1+tanh(v_tmp/10)).*(1-q_tmp)-q_tmp/tau_d_q;
    s_inc=q_tmp.*(1-s_tmp)/tau_r-s_tmp/tau_d;
    
    v_old=v;
    v=v+dt*v_inc;
    m=m_inf(v);
    h=h+dt*h_inc;
    n=n+dt*n_inc;
    q=q+dt*q_inc;
    s=s+dt*s_inc;

    % Find which neurons just fired:
    ind=find(v_old>=-20&v<-20);
    
    % Count how many there are: 
    L=length(ind);
    
    % Find the exact times at which they fired by 
    % linear interpolation between (k-1)*dt and k*dt:
    if L>0,
        t_spikes(num_spikes+1:num_spikes+L)= ...
        ((k-1)*dt*(-v(ind)-20)+k*dt*(20+v_old(ind)))...
        ./(v_old(ind)-v(ind));
        i_spikes(num_spikes+1:num_spikes+L)=ind;
        num_spikes=num_spikes+L;
    end;
end;

ind=find(i_spikes==1);
L=length(ind);

% Compute the frequency of neuron 1 at the end of the 
% simulation, and print it out for information:
frequency=1000/(t_spikes(ind(L))-t_spikes(ind(L-1)))

t=[0:m_steps]*dt;
subplot(211);
set(gca,'Fontsize',16);

% Plot spike rastergram:
plot(t_spikes,i_spikes,'.r','Markersize',10);
axis([t_final-200,t_final,0,N+1]);
shg;
hold off;
xlabel('$t$ [ms]','Fontsize',20); ylabel('neuron \#','Fontsize',20);
