function rtm_init=rtm_init(i_ext,phi_vec);



% input: i_ext=external drive (a single number)
%        phi_vec = column vector of phases at which 
%                  neurons are to be initialized

%        The length, num, of i_ext is the total number 
%        of neurons.

% output: a num-by-3 array called rtm_init. The columns
%         contain values of v, h, and n. If i_ext is below
%         the firing threshold, then the points (v,h,n)
%         (the rows of the matrix rtm_init) are all equal 
%         to each other, and equal to the stable equilibrium 
%         point. 



%         If i_ext is above the firing threshold, then the
%         (v,h,n) lie on the limit cycle, with the phases
%         given by phi_vec.

%         The function also computes the period, T. It is
%         passed back to the program using this function 
%         through a "global T". 


global T; 

num=length(phi_vec);

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

t_final=5000;           % if fewer than 5 spikes occur
                        % by this time, the program gives
                        % up and sets the (v,h,n) equal 
                        % to the values at time t_final.
                        
dt=0.001; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); v=z; m=z; h=z; n=z;

v(1)=-70;
m(1)=m_inf(v(1));
h(1)=h_inf(v(1));
n(1)=n_inf(v(1));
t=0;
k=1;

num_spikes=0;

while num_spikes<5 & t<t_final,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+ ...
        g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_l*(v_l-v(k))+i_ext)/c;
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

    v_tmp=v(k)+dt05*v_inc;   
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
 
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+ ...
           g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    
    if v(k)>=-20 & v(k+1)<-20,
        num_spikes=num_spikes+1;
        t_spikes(num_spikes)= ...
            ((k-1)*dt*(-20-v(k+1))+k*dt*(v(k)+20))/ ...
            (v(k)-v(k+1));
    end;
    t=k*dt;
    k=k+1;
end;
if num_spikes<5,
    rtm_init(:,1)=v(k);
    rtm_init(:,2)=h(k);
    rtm_init(:,3)=n(k);
    T=inf;
end;
if num_spikes==5,
    T=t_spikes(5)-t_spikes(4);
    for i=1:num,
        phi0=phi_vec(i);
        t0=phi0*T+t_spikes(4);
        k=floor(t0/dt)+1;
        rtm_init(i,1)= ...
            (v(k+1)*(t0-(k-1)*dt)+v(k)*(k*dt-t0))/dt;
        rtm_init(i,2)= ...
            (h(k+1)*(t0-(k-1)*dt)+h(k)*(k*dt-t0))/dt;
        rtm_init(i,3)= ...
            (n(k+1)*(t0-(k-1)*dt)+n(k)*(k*dt-t0))/dt;
    end;
end;
        