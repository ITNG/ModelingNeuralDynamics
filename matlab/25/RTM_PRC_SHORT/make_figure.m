clear; clf;
tic

global T; 

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

i_ext=0.30; 


N=200;
delta_v=4;

dt=0.001;
dt05=dt/2;

phi_vec=[1:N]'/N-1/(2*N);


initial_vector=rtm_init(i_ext,phi_vec);
                                % initialization in
                                % splay state
                                
v=initial_vector(:,1)+delta_v;  % At time 0, the j-th neuron is at
                                % phase (j-1/2)/N, j=1,2,...,N. 
                                % Then v is raised by delta_v, modeling
                                % a very short input pulse. We simulate
                                % until each of the neurons has fired, 
                                % then calculate by how much the firing has been
                                % advanced or delayed by the synaptic input.
m=m_inf(v);
h=initial_vector(:,2);
n=initial_vector(:,3);



% Record for each neuron how many spikes have occurred
% so far: 
num_spikes=zeros(N,1);  

% Record for each neuron at which time those spikes occurred:
t_spikes=[]; 


k=0;    % number of time steps done so far


while min(num_spikes)<1,
    
    k=k+1;
    
    v_inc=(g_k*n.^4.*(v_k-v)+g_na*m.^3.*h.*(v_na-v)+...
        g_l*(v_l-v)+i_ext)/c;
    h_inc=alpha_h(v).*(1-h)-beta_h(v).*h;
    n_inc=alpha_n(v).*(1-n)-beta_n(v).*n;
    
    
    v_tmp=v+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h+dt05*h_inc;
    n_tmp=n+dt05*n_inc;
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
        g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c; 
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v_old=v;
    v=v+dt*v_inc;
    m=m_inf(v);
    h=h+dt*h_inc;
    n=n+dt*n_inc;

    ind=find(v_old>=-20&v<-20); % Determine which neurons
                                % fired in this time step,
                                % if any.
    L=length(ind);              % L is the number of neurons
                                % that just fired. 
    for i=1:L,
        num_spikes(ind(i))=num_spikes(ind(i))+1;
        t_spikes(ind(i),num_spikes(ind(i)))= ...
        ((k-1)*dt*(-20-v(ind(i)))+k*dt*(20+v_old(ind(i))))...
        /(-v(ind(i))+v_old(ind(i))); 
    end;

end;

t_star=t_spikes(1:N,1);

g_vec=-t_star/T-phi_vec+1;


subplot(111);
plot(phi_vec,g_vec,'-k','Linewidth',6);
set(gca,'Fontsize',24);
axis([0,1,0,1]);
axis('square');
xlabel('$\varphi$','Fontsize',32);
ylabel('$g$','Fontsize',32);
hold on;
plot([0,1],[1,0],'--k','Linewidth',2);
hold off;
shg;


