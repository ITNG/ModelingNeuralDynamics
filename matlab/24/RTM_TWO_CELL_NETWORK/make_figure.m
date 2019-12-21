clear; clf;

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

i_ext=0.3; 

tau_r=0.5; tau_peak=0.5; tau_d=2;
tau_d_q=tau_d_q_function(tau_d,tau_r,tau_peak);


N=2;
g_syn=0.0075*29;
t_final=1000;

dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

% Initialize the two neurons at the same phase, arbitrarily
% taken to be 0.4 here:
initial_vector=rtm_init(i_ext,[0.4;0.4]);

v=initial_vector(:,1);
m=m_inf(v);
h=initial_vector(:,2);
n=initial_vector(:,3);
q=zeros(N,1);
s=zeros(N,1);

num_spikes=0;   % number of spikes so far
t_spikes=[];    % times of those spikes
i_spikes=[];    % The i-th spike is a spike of 
                % neuron number i_spikes(i). 

done=0;         % "done" will be 0 until the 
                % small perturbation from synchrony
                % has been introduced, and it will then
                % be set to 1. 


for k=1:m_steps,
    
    v_inc=(g_k*n.^4.*(v_k-v)+g_na*m.^3.*h.*(v_na-v)+...
        g_l*(v_l-v)-g_syn*(sum(s)-s).*v+i_ext)/c;
    h_inc=alpha_h(v).*(1-h)-beta_h(v).*h;
    n_inc=alpha_n(v).*(1-n)-beta_n(v).*n;
    q_inc=5*(1+tanh(v/10)).*(1-q)-q/tau_d;
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
    
    % Find out which neurons just fired:
    ind=find(v_old>=-20&v<-20);
    
    % Count how many such neurons there are:
    L=length(ind);
    
    % Compute exactly when, in the interval 
    % from (k-1)*dt to k*dt, they fired, by
    % linear interpolation:
    if L>0,
        t_spikes(num_spikes+1:num_spikes+L)= ...
        ((k-1)*dt*(-v(ind)-20)+k*dt*(20+v_old(ind)))...
        ./(-v(ind)+v_old(ind));
        i_spikes(num_spikes+1:num_spikes+L)=ind;
        num_spikes=num_spikes+L;
    end;
    
    % If the 10-th action potential (5 for each neuron)
    % has happened, and if the perturbation has not yet 
    % been introduced, introduce it now:
    if num_spikes==10 & done==0,
        v(1)=v(1)-10^(-5);
        done=1;
    end;
end;

ind=find(i_spikes==1);
t1=t_spikes(ind);
ind=find(i_spikes==2);
t2=t_spikes(ind);
subplot(211);
plot(t2-t1,'.k','Markersize',20);
set(gca,'Fontsize',16);
xlabel('spike \#', 'Fontsize',20);
title('spike time difference [ms]','Fontsize',20);
axis([0,25,-2,2]);

hold on;
h=plot([5,5],[-3,-4]);

% This command allows you to plot outside the plotting
% window (when you look at the plot that this code
% generates, you'll see why I wanted that):
set(h,'clipping','off');

h=text(0,-4.5,'time when $v_1$ is lowered','Fontsize',20);
arrow(0,25,-2,6,5,-3,[0;1],0.025,2,'-k');
hold off;

shg;

