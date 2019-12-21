clear; 

tic 
format long;
    
c=1;
g_k=80; 
g_na=100;
g_l=0.1; 
v_k=-100;
v_na=50;
v_l=-67;

epsilon=10^(-5);
                     

dt=0.01; dt05=dt/2;

z=zeros(round(5000/dt)+1,1);    % allocate space for v, m, h, n
v=z; m=z; h=z; n=z;

w=zeros(5,1);

for ijk=0:4,
    g_bar=0.15+ijk*0.05  % strength of tonic inhibition
    v_rev=-75;         % reversal potential of tonic inhibition;
    g_mean=mean(g((1:10000)/10000*25));
    g_factor=g_bar/g_mean; % factor by which g(t) must be multiplied
                           % to make the mean equal to g_bar. 

    v(1)=-70; 
    m(1)=m_inf(v(1));
    h(1)=0.7; 
    n(1)=0.6; 

    done=0;
    t_final=5000;
    m_steps=round(t_final/dt);
    g_store=g([0:m_steps]*dt)*g_factor;

    i_ext_left=0; i_ext_right=3;

    while (i_ext_right-i_ext_left)/((i_ext_right+i_ext_left)/2)>epsilon,
        i_ext=(i_ext_left+i_ext_right)/2;
        num_spikes=0;
        for k=1:m_steps,

            t_old=(k-1)*dt;
            t_tmp=(k-1/2)*dt;

            v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
                g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
                g_store(k)*(v_rev-v(k))+ ...
                i_ext)/c;
            h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
            n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

            v_tmp=v(k)+dt05*v_inc;
            m_tmp=m_inf(v_tmp);
            h_tmp=h(k)+dt05*h_inc;
            n_tmp=n(k)+dt05*n_inc;

            v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
                g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
                (g_store(k)+g_store(k+1))*0.5*(v_rev-v_tmp)+ ...
                i_ext)/c;
            h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
            n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;

            v(k+1)=v(k)+dt*v_inc;
            m(k+1)=m_inf(v(k+1));
            h(k+1)=h(k)+dt*h_inc;
            n(k+1)=n(k)+dt*n_inc;

            if v(k+1)<-20 && v(k)>=-20,
                num_spikes=num_spikes+1;
                t_spikes(num_spikes)= ...
                    (k*dt*(20+v(k))+(k-1)*dt*(-20-v(k+1)))/(v(k)-v(k+1));
            end;

        end;
        f=num_spikes/t_final*1000;
        if f==0, 
            i_ext_left=i_ext;
        else
            i_ext_right=i_ext;
        end;
    end;
    I_L=(i_ext_left+i_ext_right)/2

    i_ext_left=0; i_ext_right=3;

    while (i_ext_right-i_ext_left)/((i_ext_right+i_ext_left)/2)>epsilon,
        i_ext=(i_ext_left+i_ext_right)/2;
        num_spikes=0;
        for k=1:m_steps,

            t_old=(k-1)*dt;
            t_tmp=(k-1/2)*dt;

            v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
                g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
                g_store(k)*(v_rev-v(k))+ ...
                i_ext)/c;
            h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
            n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

            v_tmp=v(k)+dt05*v_inc;
            m_tmp=m_inf(v_tmp);
            h_tmp=h(k)+dt05*h_inc;
            n_tmp=n(k)+dt05*n_inc;

            v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
                g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
                (g_store(k)+g_store(k+1))*0.5*(v_rev-v_tmp)+ ...
                i_ext)/c;
            h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
            n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;

            v(k+1)=v(k)+dt*v_inc;
            m(k+1)=m_inf(v(k+1));
            h(k+1)=h(k)+dt*h_inc;
            n(k+1)=n(k)+dt*n_inc;

            if v(k+1)<-20 && v(k)>=-20,
                num_spikes=num_spikes+1;
                t_spikes(num_spikes)= ...
                    (k*dt*(20+v(k))+(k-1)*dt*(-20-v(k+1)))/(v(k)-v(k+1));
            end;

        end;
        f=num_spikes/t_final*1000;
        if f<39, 
            i_ext_left=i_ext;
        else
            i_ext_right=i_ext;
        end;
    end;
    I_R=(i_ext_left+i_ext_right)/2
    w(ijk+1)=I_R-I_L
end;
w

ratios=w(1:4)./w(2:5)

toc

