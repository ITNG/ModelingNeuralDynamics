clear; clf;

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

i_ext=0.30; 

tau_r=0.5; tau_peak=0.5; tau_d=2;
tau_d_q=tau_d_q_function(tau_d,tau_r,tau_peak);


N=30;
g_syn=0.0075;
t_final=10000;


dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

initial_vector=rtm_init(i_ext,[N-1:-1:0]'/N+1/(2*N));
                                % initialization in
                                % splay state
v=initial_vector(:,1);
m=m_inf(v);
h=initial_vector(:,2);
n=initial_vector(:,3);
q=zeros(N,1);
s=zeros(N,1);

num_spikes=0;       % counts the number of spikes
                    % that have occurred
t_spikes=[];        % records spike times
i_spikes=[];        % The i-th spike is a spike of neuron
                    % i_spikes(i). (It occurs at time
                    % t_spikes(i).) 

for k=1:m_steps,
    
    v_inc=(g_k*n.^4.*(v_k-v)+g_na*m.^3.*h.*(v_na-v)+...
        g_l*(v_l-v)-g_syn*(sum(s)-s).*v+i_ext)/c;
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
        g_syn*(sum(s_tmp)*v_tmp-s_tmp.*v_tmp)+i_ext)/c;
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

    ind=find(v_old>=-20&v<-20); % Determine which neurons
                                % fired in this time step,
                                % if any.
    L=length(ind);              % L is the number of neurons
                                % that just fired. 
    if L>0,
        t_spikes(num_spikes+1:num_spikes+L)= ...
        ((k-1)*dt*(-v(ind)-20)+k*dt*(20+v_old(ind)))...
        ./(-v(ind)+v_old(ind)); % Find the time at which 
                                % v crossed -20. This time
                                % lies between (k-1)*dt and
                                % k*dt, and is found here
                                % by linear interpolation. 
        i_spikes(num_spikes+1:num_spikes+L)=ind;
        num_spikes=num_spikes+L;
    end;
end;

ind=find(i_spikes==1);
L=length(ind);
frequency=1000/(t_spikes(ind(L))-t_spikes(ind(L-1)))
                                % This is the frequency 
                                % of the first neuron, 
                                % computed based on the 
                                % last interspike interval,
                                % and printed out for 
                                % information.
t=[0:m_steps]*dt;
subplot(211);
set(gca,'Fontsize',16);
plot(t_spikes,i_spikes,'.r','Markersize',10);
axis([t_final-200,t_final,0,N+1]);
shg;
hold off;
xlabel('$t$ [ms]','Fontsize',20); 
ylabel('neuron \#','Fontsize',20);
    
    
