function olm_init=olm_init(i_ext,phi_vec)



% input: i_ext=column vector of external drives 
%        phi_vec = column vector of phases at which 
%                  neurons are to be initialized

%        The length, num, of i_ext is the total number 
%        of neurons.

% output: a num-by-6 array called olm_init. The columns
%         contain values of v, h, n, r, a, b. If i_ext(i) is below
%         the firing threshold, then the i-th row of the
%         matrix olm_init contains the stable equilibrium point.
%         If i_ext(i) is above the firing threshold, then the
%         i-th row of olm_init is a point (v,h,n,r,a,b) on the limit
%         cycle, at phase phi_vec(i). 


num=length(phi_vec);
olm_init=zeros(num,6);

c=1.3;
g_k=23; 
g_na=30;
g_l=0.05;
v_k=-100;
v_na=90;
v_l=-70;
g_h=12;
g_A=22;
v_h=-32.9;
v_A=-90;

max_spikes=3;
t_final=2000;           % if fewer than max_spikes spikes occur
                        % by this time, the program gives
                        % up and sets (v,h,n,r,a,b) equal 
                        % to the values at time t_final.
                        
dt=0.01; dt05=dt/2; 

v=-70*ones(num,1);
m=m_o_inf(v);
h=h_o_inf(v);
n=n_o_inf(v);
r=r_o_inf(v);
a=a_o_inf(v);
b=b_o_inf(v);

t=0;

num_spikes=zeros(num,1); 
                                       
done=zeros(num,1);                     % done(i)=1 indicates that we are 
                                       % done with neuron i.
                                       

t_spikes=zeros(num,max_spikes);

while sum(done)<num && t<t_final,
    
    v_old=v;
    h_old=h;
    n_old=n;
    r_old=r;
    a_old=a;
    b_old=b;
    t_old=t;
    
    v_inc=(g_k*n.^4.*(v_k-v)+ ...
        g_na*m.^3.*h.*(v_na-v)+ ...
        g_h*r.*(v_h-v)+ ...
        g_A*a.*b.*(v_A-v)+ ...
        g_l*(v_l-v)+i_ext)/c;
    h_inc=(h_o_inf(v)-h)./tau_h_o(v);
    n_inc=(n_o_inf(v)-n)./tau_n_o(v);
    r_inc=(r_o_inf(v)-r)./tau_r_o(v);
    a_inc=(a_o_inf(v)-a)./tau_a_o(v);
    b_inc=(b_o_inf(v)-b)./tau_b_o(v);
   
    v_tmp=v+dt05*v_inc;
    m_tmp=m_o_inf(v);
    h_tmp=h+dt05*h_inc;
    n_tmp=n+dt05*n_inc;
    r_tmp=r+dt05*r_inc;
    a_tmp=a+dt05*a_inc;
    b_tmp=b+dt05*b_inc;
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+ ...
                 g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
                 g_l*(v_l-v_tmp)+ ...
                 g_h*r_tmp.*(v_h-v_tmp)+ ...
                 g_A*a_tmp.*b_tmp.*(v_A-v_tmp)+ ...
                 i_ext)/c;
    h_inc=(h_o_inf(v_tmp)-h_tmp)./tau_h_o(v_tmp);
    n_inc=(n_o_inf(v_tmp)-n_tmp)./tau_n_o(v_tmp);
    r_inc=(r_o_inf(v_tmp)-r_tmp)./tau_r_o(v_tmp);
    a_inc=(a_o_inf(v_tmp)-a_tmp)./tau_a_o(v_tmp);
    b_inc=(b_o_inf(v_tmp)-b_tmp)./tau_b_o(v_tmp);
    
    v=v+dt*v_inc;
    m=m_o_inf(v);
    h=h+dt*h_inc;
    n=n+dt*n_inc;
    r=r+dt*r_inc;
    a=a+dt*a_inc;
    b=b+dt*b_inc;
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
        olm_init(k,1)=(v_old(k)*(t-thr(k))+v(k)*(thr(k)-t_old))/dt;
        olm_init(k,2)=(h_old(k)*(t-thr(k))+h(k)*(thr(k)-t_old))/dt;
        olm_init(k,3)=(n_old(k)*(t-thr(k))+n(k)*(thr(k)-t_old))/dt;
        olm_init(k,4)=(r_old(k)*(t-thr(k))+r(k)*(thr(k)-t_old))/dt;
        olm_init(k,5)=(a_old(k)*(t-thr(k))+a(k)*(thr(k)-t_old))/dt;
        olm_init(k,6)=(b_old(k)*(t-thr(k))+b(k)*(thr(k)-t_old))/dt;
    end;
    done(ind)=1;
    
end;
ind=find(done==0);
olm_init(ind,1)=v(ind);
olm_init(ind,2)=h(ind);
olm_init(ind,3)=n(ind);
olm_init(ind,4)=r(ind);
olm_init(ind,5)=a(ind);
olm_init(ind,6)=b(ind);







        