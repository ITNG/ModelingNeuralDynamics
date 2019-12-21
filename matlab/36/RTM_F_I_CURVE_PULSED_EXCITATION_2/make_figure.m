clear; clf;

global alpha Period ave;        % Parameters determining the shape and
                                % period of the input pulses. See below
                                % and also shape.m. 

alpha=1;        % determines sharpness of input pulses
Period=25;      % period of input pulses
N=2000; 
ave=mean(exp(alpha*cos(pi*[0:N-1]/N).^2)-1);    % see shape.m



tic 

c=1;
g_k=80; 
g_na=100;
g_l=0.2; 
v_k=-100;
v_na=50;
v_l=-67;


i_ext_low=0; i_ext_high=2; 
i_ext_vec=i_ext_low+(0:200)/200*(i_ext_high-i_ext_low); % external drives
                                                        % to be considered


dt=0.01; dt05=dt/2;
t_final=1000;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);    % allocate space for v, m, h, n
v=z; m=z; h=z; n=z;
t=(0:m_steps)*dt;
shape_store=shape(t);              % The input pulses always have the 
                                   % same shape, but with varying 
                                   % amplitude. shape_store has temporal
                                   % average 1. 

i_ext=i_ext_vec(1);

v(1)=-70; 
m(1)=m_inf(v(1));
h(1)=0.7; 
n(1)=0.6; 

num_spikes=0;       % number of spikes simulated so far                
k=1;                % in the arrays v, m, h, n, entry k is the latest
                    % one that has been computed. 
done=0;             % done=0 indicates that for the current value
                    % of i_ext, we aren't done yet. 

while done==0,
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
        i_ext)/c;
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
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
    
    if num_spikes==2,                       % if 2 spikes have occurred,
        f=1000/(t_spikes(2)-t_spikes(1));   % take the time between them
        done=1;                             % to be the period. 
    end;
    k=k+1;
    if k==m_steps+1,                        % if time t_final has been
        f=0;                                % reached, take the period
        done=1;                             % to be zero. (If two spikes
    end;                                    % had occurred prior to time
end;                                        % t_final, we would never
                                            % have reached time t_final.)
f_vec(1)=f;
ijk=1

for ijk=2:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);

    v(1)=v(k); 
    m(1)=m(k);
    h(1)=h(k);
    n(1)=n(k);

    done=0;
    num_spikes=0;
    k=1;
    
    while done==0,
    
        v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
            g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
            i_ext)/c;
        h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
        n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

        v_tmp=v(k)+dt05*v_inc;
        m_tmp=m_inf(v_tmp);
        h_tmp=h(k)+dt05*h_inc;
        n_tmp=n(k)+dt05*n_inc;

        v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
            g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
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

        if num_spikes==2,
            f=1000/(t_spikes(2)-t_spikes(1));
            done=1;
        end;
        k=k+1;
        if k==m_steps+1,
            f=0;
            done=1;
        end;
        
    end;
    f_vec(ijk)=f;
    ijk
end;

subplot(211);
plot(i_ext_vec,f_vec,'.r','Markersize',10);
set(gca,'Fontsize',16);
xlabel('$I$ [$\mu$A/cm$^2$]','Fontsize',20);
ylabel('$f$','Fontsize',20);
hold on;


i_ext=i_ext_vec(1);
i_ext_p=i_ext*shape_store;      % this is the periodic excitatory input
                                % with time average i_ext.
v(1)=-70; 
m(1)=m_inf(v(1));
h(1)=0.7; 
n(1)=0.6; 


num_spikes=0;

for k=1:m_steps,
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
        i_ext_p(k))/c;
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
        (i_ext_p(k)+i_ext_p(k+1))/2)/c;
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

f_vec(1)=num_spikes/t_final*1000;
ijk=1

for ijk=2:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    i_ext_p=i_ext*shape_store;

    v(1)=v(k); 
    m(1)=m(k);
    h(1)=h(k);
    n(1)=n(k);

    num_spikes=0;

    for k=1:m_steps,
    
        v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
            g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
            i_ext_p(k))/c;
        h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
        n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

        v_tmp=v(k)+dt05*v_inc;
        m_tmp=m_inf(v_tmp);
        h_tmp=h(k)+dt05*h_inc;
        n_tmp=n(k)+dt05*n_inc;

        v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
            g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
            (i_ext_p(k)+i_ext_p(k+1))/2)/c;
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
    f_vec(ijk)=num_spikes/t_final*1000;
    ijk
end;

plot(i_ext_vec,f_vec,'.b','Markersize',10);
hold off;


shg;

toc

