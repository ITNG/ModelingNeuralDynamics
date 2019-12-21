clear; clf;

c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;

i_ext_vec=3+[0:60]/60*10;   % vector of external drives for which calculations
                            % will be done.

dt=0.01; dt05=dt/2;

z=zeros(round(5000/dt),1);    
v=z; m=z; h=z; n=z;         % pre-allocate space for v, m, h, n, for
                            % 5 seconds of simulated time.

i_ext=i_ext_vec(1);

v(1)=-70; 
m(1)=alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
h(1)=0.7; 
n(1)=0.6; 

done=0;                 % done=0 as long as we aren't done. 
num_spikes=0;           % number of action potentials so far.
N=round(1000/dt);       % number of time steps in a second.
k=1;                    % k-1 time steps have been done.

while done==0,
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    m_inc=alpha_m(v(k))*(1-m(k))-beta_m(v(k))*m(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m(k)+dt05*m_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m(k)+dt*m_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    
    if mod(k-1,N)==0 & k>1,
        maxv=max(v(k-N+1:k+1));
        minv=min(v(k-N+1:k+1));
        maxm=max(m(k-N+1:k+1));
        minm=min(m(k-N+1:k+1));
        maxh=max(h(k-N+1:k+1));
        minh=min(h(k-N+1:k+1));
        maxn=max(n(k-N+1:k+1));
        minn=min(n(k-N+1:k+1));
        if (maxv-minv)<0.0001*abs(maxv+minv) & ...
           (maxm-minm)<0.0001*abs(maxm+minm) & ...
           (maxh-minh)<0.0001*abs(maxh+minh) & ...
           (maxn-minn)<0.0001*abs(maxn+minn),
            f=0;
            done=1;
        end;
    end;
    
    if v(k+1)<-20 & v(k)>=-20,
        num_spikes=num_spikes+1;
        t_spikes(num_spikes)= ...
            (k*dt*(v(k)+20)+(k-1)*dt*(-20-v(k+1)))/(v(k)-v(k+1));
    end;
    
    if num_spikes==4,
        f=1000/(t_spikes(4)-t_spikes(3));
        done=1;
    end;
    k=k+1;
end;


f_vec(1)=f;

for ijk=2:length(i_ext_vec),    % upward sweep
    i_ext=i_ext_vec(ijk);

    v(1)=v(k);          % set phase space position in the location where
    m(1)=m(k);          % previous simulation ended.
    h(1)=h(k);
    n(1)=n(k);

    done=0;
    num_spikes=0;
    k=1;
    while done==0,
    
        v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
            g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
        m_inc=alpha_m(v(k))*(1-m(k))-beta_m(v(k))*m(k);
        h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
        n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

        v_tmp=v(k)+dt05*v_inc;
        m_tmp=m(k)+dt05*m_inc;
        h_tmp=h(k)+dt05*h_inc;
        n_tmp=n(k)+dt05*n_inc;

        v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
            g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
        m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
        h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
        n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;

        v(k+1)=v(k)+dt*v_inc;
        m(k+1)=m(k)+dt*m_inc;
        h(k+1)=h(k)+dt*h_inc;
        n(k+1)=n(k)+dt*n_inc;

        if mod(k-1,N)==0 & k>1,
            maxv=max(v(k-N+1:k+1));
            minv=min(v(k-N+1:k+1));
            maxm=max(m(k-N+1:k+1));
            minm=min(m(k-N+1:k+1));
            maxh=max(h(k-N+1:k+1));
            minh=min(h(k-N+1:k+1));
            maxn=max(n(k-N+1:k+1));
            minn=min(n(k-N+1:k+1));
            if (maxv-minv)<0.0001*abs(maxv+minv) & ...
               (maxm-minm)<0.0001*abs(maxm+minm) & ...
               (maxh-minh)<0.0001*abs(maxh+minh) & ...
               (maxn-minn)<0.0001*abs(maxn+minn),
                f=0;
                done=1;
            end;
        end;

        if v(k+1)<-20 & v(k)>=-20,
            num_spikes=num_spikes+1;
            t_spikes(num_spikes)= ...
                (k*dt*(v(k)+20)+(k-1)*dt*(-20-v(k+1)))/(v(k)-v(k+1));
        end;

        if num_spikes==4,
            f=1000/(t_spikes(4)-t_spikes(3));
            done=1;
        end;
        k=k+1;
    end;
    f_vec(ijk)=f;   
    ijk         % tell us how many external drives have been handled yet, 
                % so we don't get too impatient while the code is running.
end;    

subplot(211);
plot(i_ext_vec(1:2:length(i_ext_vec)),f_vec(1:2:length(f_vec)), ...
    '.k','Markersize',15);          % plot only 31 of the 61 computed
                                    % data points, for neatness


ind=find(f_vec==0);
ind=max(ind);
I_c=(i_ext_vec(ind)+i_ext_vec(ind+1))/2;
hold on;
plot([i_ext_vec(ind+1),i_ext_vec(ind+1)],[0,f_vec(ind+1)], ...
    '--b','Linewidth',1);


f_vec(length(i_ext_vec))=f;
for ijk=length(i_ext_vec)-1:-1:1,   % downward sweep
    i_ext=i_ext_vec(ijk);

    v(1)=v(k);  % start where the previous simulation ended
    m(1)=m(k);
    h(1)=h(k);
    n(1)=n(k);

    done=0;
    num_spikes=0;
    k=1;
    while done==0,
    
        v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
            g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
        m_inc=alpha_m(v(k))*(1-m(k))-beta_m(v(k))*m(k);
        h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
        n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);

        v_tmp=v(k)+dt05*v_inc;
        m_tmp=m(k)+dt05*m_inc;
        h_tmp=h(k)+dt05*h_inc;
        n_tmp=n(k)+dt05*n_inc;

        v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
            g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
        m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
        h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
        n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;

        v(k+1)=v(k)+dt*v_inc;
        m(k+1)=m(k)+dt*m_inc;
        h(k+1)=h(k)+dt*h_inc;
        n(k+1)=n(k)+dt*n_inc;

        if mod(k-1,N)==0 & k>1,
            maxv=max(v(k-N+1:k+1));
            minv=min(v(k-N+1:k+1));
            maxm=max(m(k-N+1:k+1));
            minm=min(m(k-N+1:k+1));
            maxh=max(h(k-N+1:k+1));
            minh=min(h(k-N+1:k+1));
            maxn=max(n(k-N+1:k+1));
            minn=min(n(k-N+1:k+1));
            if (maxv-minv)<0.0001*abs(maxv+minv) & ...
               (maxm-minm)<0.0001*abs(maxm+minm) & ...
               (maxh-minh)<0.0001*abs(maxh+minh) & ...
               (maxn-minn)<0.0001*abs(maxn+minn),
                f=0;
                done=1;
            end;
        end;

        if v(k+1)<-20 & v(k)>=-20,
            num_spikes=num_spikes+1;
            t_spikes(num_spikes)= ...
                (k*dt*(v(k)+20)+(k-1)*dt*(-20-v(k+1)))/(v(k)-v(k+1));
        end;

        if num_spikes==4,
            f=1000/(t_spikes(4)-t_spikes(3));
            done=1;
        end;
        k=k+1;
    end;
    f_vec(ijk)=f;
    ijk         % tell us which external drive we are currently testing,
                % so we don't get too impatient while the code is running.
end;
plot(i_ext_vec(1:2:length(i_ext_vec)),f_vec(1:2:length(f_vec)), ...
    'ok','Markersize',10,'Linewidth',1); % plot only 31 of the 61 computed
                                         % data points, for neatness

ind=find(f_vec==0);
ind=max(ind);
I_star=(i_ext_vec(ind)+i_ext_vec(ind+1))/2;

I_c
I_star

set(gca,'Xtick',[5,8,11]);
plot([i_ext_vec(ind+1),i_ext_vec(ind+1)],[0,f_vec(ind+1)], ...
    '--b','Linewidth',1);

axis([3,13,0,100]);
hold off;
set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$f$','Fontsize',20);
text(I_star-0.1,-15,'$I_\ast$','Fontsize',20,'Color','b');
text(I_c-0.1,-15,'$I_c$','Fontsize',20,'Color','b');
shg;