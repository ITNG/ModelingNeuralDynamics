function rtm_init_with_m_current=rtm_init_with_m_current(i_ext,phi_vec,g_m)



% input: i_ext=column vector of external drives 
%        phi_vec = column vector of phases at which 
%                  neurons are to be initialized
%        g_m = strength of m-current

%        The length, num, of i_ext is the total number 
%        of neurons.

% output: a num-by-4 array called rtm_init_with_m_current. The columns
%         contain values of v, h, n, and w. If i_ext(i) is below
%         the firing threshold, then the i-th row of the
%         matrix rtm_init_with_m_current contains the stable equilibrium point.


%         If i_ext(i) is above the firing threshold, then the
%         i-th row of rtm_init_with_m_current is a point (v,h,n,w) on the limit
%         cycle, at phase phi_vec(i). 



num=length(phi_vec);
rtm_init_with_m_current=zeros(num,4);

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

max_spikes=3;
t_final=2000;           % if fewer than max_spikes spikes occur
                        % by this time, the program gives
                        % up and sets (v,h,n) equal 
                        % to the values at time t_final.
                        
dt=0.01; dt05=dt/2; 

v=-70*ones(num,1);
m=m_e_inf(v);
h=h_e_inf(v);
n=n_e_inf(v);
w=zeros(num,1);

t=0;

num_spikes=zeros(num,1); 
                                       
done=zeros(num,1);                     % done(i)=1 indicates that we are 
                                       % done with neuron i.
                                       

t_spikes=zeros(num,max_spikes);

while sum(done)<num && t<t_final,
    
    v_old=v;
    h_old=h;
    n_old=n;
    w_old=w;
    t_old=t;
    
    v_inc=(g_k*n.^4.*(v_k-v)+ ...
        g_na*m.^3.*h.*(v_na-v)+ ...
        g_l*(v_l-v)+g_m*w.*(v_k-v)+i_ext)/c;
    h_inc=(h_e_inf(v)-h)./tau_h_e(v);
    n_inc=(n_e_inf(v)-n)./tau_n_e(v);
    w_inc=(w_inf(v)-w)./tau_w(v);
   
    v_tmp=v+dt05*v_inc;
    m_tmp=m_e_inf(v);
    h_tmp=h+dt05*h_inc;
    n_tmp=n+dt05*n_inc;
    w_tmp=w+dt05*w_inc;
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
                 g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
                 g_l*(v_l-v_tmp)+g_m*w_tmp.*(v_k-v_tmp)+i_ext)/c;
    h_inc=(h_e_inf(v_tmp)-h_tmp)./tau_h_e(v_tmp);
    n_inc=(n_e_inf(v_tmp)-n_tmp)./tau_n_e(v_tmp);
    w_inc=(w_inf(v_tmp)-w_tmp)./tau_w(v_tmp);
    
    v=v+dt*v_inc;
    m=m_e_inf(v);
    h=h+dt*h_inc;
    n=n+dt*n_inc;
    w=w+dt*w_inc;
    t=t+dt;
    
    ind=find(v_old>=-20 & v<-20);
    l=length(ind);
    for i=1:l,
        k=ind(i);
        num_spikes(k)=num_spikes(k)+1;
        t_spikes(k,num_spikes(k))= ...
            (t_old*(-20-v(k))+t*(v_old(k)-(-20)))/ ...
            (v_old(k)-v(k));
    end;
    
    thr=t_spikes(:,max_spikes)+ ...
        phi_vec.*(t_spikes(:,max_spikes)-t_spikes(:,max_spikes-1));
    ind=find(num_spikes==max_spikes & t>thr & t_old <= thr);
    l=length(ind);
    
    for i=1:l,
        k=ind(i);
        rtm_init_with_m_current(k,1)= ...
            (v_old(k)*(t-thr(k))+v(k)*(thr(k)-t_old))/dt;
        rtm_init_with_m_current(k,2)= ...
            (h_old(k)*(t-thr(k))+h(k)*(thr(k)-t_old))/dt;
        rtm_init_with_m_current(k,3)= ...
            (n_old(k)*(t-thr(k))+n(k)*(thr(k)-t_old))/dt;
        rtm_init_with_m_current(k,4)= ...
            (w_old(k)*(t-thr(k))+w(k)*(thr(k)-t_old))/dt;
    end;
    done(ind)=1;
    
end;
ind=find(done==0);
rtm_init_with_m_current(ind,1)=v(ind);
rtm_init_with_m_current(ind,2)=h(ind);
rtm_init_with_m_current(ind,3)=n(ind);
rtm_init_with_m_current(ind,4)=w(ind);







        