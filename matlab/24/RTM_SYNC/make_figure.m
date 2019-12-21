clear; clf;
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

i_ext=0.3; 

N=30;           % number of neurons in the network

t_final=200;
dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

v=-70*ones(N,1);
m=m_inf(v);
h=h_inf(v);
n=n_inf(v);


num_spikes=0;   % number of spikes so far
t_spikes=[];    % times of those spikes
i_spikes=[];    % The i-th spike is a spike of neuron
                % number i_spikes(i). 

for k=1:m_steps,
    
    v_inc=(g_k*n.^4.*(v_k-v)+g_na*m.^3.*h.*(v_na-v)+...
        g_l*(v_l-v)+i_ext)/c;
    n_inc=alpha_n(v).*(1-n)-beta_n(v).*n;
    h_inc=alpha_h(v).*(1-h)-beta_h(v).*h;
    
    v_tmp=v+dt05*v_inc;
    h_tmp=h+dt05*h_inc;
    n_tmp=n+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
        g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v_old=v;
    v=v+dt*v_inc;
    h=h+dt*h_inc;
    n=n+dt*n_inc;
    m=m_inf(v);
    
    % Find which neurons just crossed v=-20 from above:
    ind=find(v_old>=-20&v<-20);
    
    % Count how many of those neurons there are:
    L=length(ind);
    
    % Find exactly when they crossed v=-20 by linear 
    % interpolation between times (k-1)*dt and k*dt:
    if L>0,
        t_spikes(num_spikes+1:num_spikes+L)= ...
        ((k-1)*dt*(-20-v(ind))+k*dt*(v_old(ind)+20))...
        ./(v_old(ind)-v(ind));
        i_spikes(num_spikes+1:num_spikes+L)=ind;
        num_spikes=num_spikes+L;
    end;
    
end;

t=[0:m_steps]*dt;
subplot(211);
set(gca,'Fontsize',16);
plot(t_spikes,i_spikes,'.r','Markersize',10);
axis([0,t_final,0,N+1]);
shg;
hold off;
xlabel('$t$ [ms]','Fontsize',20); ylabel('neuron \#','Fontsize',20);
    
    
