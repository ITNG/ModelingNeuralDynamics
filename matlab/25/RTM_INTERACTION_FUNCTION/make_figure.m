clear; clf;
tic

global T;   % Period of model neuron. This is computed by rtm_init. 

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

N=200;
g_syn=0.1;

dt=0.01;
dt05=dt/2;

phi_vec=[1:N]'/N-1/(2*N);


initial_vector=rtm_init(i_ext,phi_vec);
                                % initialization in
                                % splay state
v=initial_vector(:,1);
m=m_inf(v);
h=initial_vector(:,2);
n=initial_vector(:,3);
q=ones(N,1);
s=zeros(N,1);


% Record for each neuron how many spikes have occurred
% so far: 
num_spikes=zeros(N,1);  

% Record for each neuron at which time those spikes occurred:
t_spikes=[]; 


k=0;    % number of time steps done so far

while min(num_spikes)<1,
    
    k=k+1;
    
    v_inc=(g_k*n.^4.*(v_k-v)+g_na*m.^3.*h.*(v_na-v)+...
        g_l*(v_l-v)- ...
        g_syn*s.*v+i_ext)/c;
    h_inc=alpha_h(v).*(1-h)-beta_h(v).*h;
    n_inc=alpha_n(v).*(1-n)-beta_n(v).*n;
    q_inc=-q/tau_d_q;
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
        g_syn*s_tmp.*v_tmp +i_ext)/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    q_inc=-q_tmp/tau_d_q;
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
    for i=1:L,
        num_spikes(ind(i))=num_spikes(ind(i))+1;
        t_spikes(ind(i),num_spikes(ind(i)))= ...
        ((k-1)*dt*(-20-v(ind(i)))+k*dt*(v_old(ind(i))+20))...
        /(v_old(ind(i))-v(ind(i))); 
    end;

end;

t_star=t_spikes(1:N,1);

g_vec=-t_star/T+1-phi_vec;


subplot(111);
plot(phi_vec,phi_vec+g_vec,'-k','Linewidth',6);
set(gca,'Fontsize',24);
axis([0,1,0,1]);
axis('square');
xlabel('$\varphi$','Fontsize',32);
ylabel('$f$','Fontsize',32);
hold on;
plot([0,1],[0,1],'--k','Linewidth',2);
hold off;
shg;
toc


