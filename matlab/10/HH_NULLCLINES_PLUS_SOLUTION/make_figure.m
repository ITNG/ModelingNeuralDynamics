clear;
v=[-100:50]+0.001;


subplot(111);
plot(v,alpha_n(v)./(alpha_n(v)+beta_n(v)),'-k','Linewidth',3);
set(gca,'Fontsize',24);
xlabel('$v$ [mV]','Fontsize',32); ylabel('$n$','Fontsize',32);
axis([-100,50,0,1]);
	
c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;
i_ext=10;
v_vec=[-100:50];

for ij=1:length(v_vec);
    v=v_vec(ij);
    n_L=0; n_R=1;
    while n_R-n_L>10^(-10),
        n=(n_L+n_R)/2;
        i_ion=g_na*m_inf(v)^3*(0.83-n)*(v_na-v)+g_k*n^4*(v_k-v)+g_l*(v_l-v)+i_ext;
        if i_ion<0,
            n_R=n;
        else
            n_L=n;
        end;
    end;
    n_vec(ij)=(n_L+n_R)/2;
    if n_vec(ij)>1-10^(-10), n_vec(ij)=2; end;
    if n_vec(ij)<10^(-10), n_vec(ij)=-1; end;
end;

hold on;
plot(v_vec,n_vec,'-r','Linewidth',3);

t_final=50;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);

v(1)=-70;
n(1)=0;
m(1)=m_inf(v(1));
h(1)=0.83-n(1);

for k=1:m_steps,
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=0.83-n_tmp;
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=0.83-n(k+1);
    if abs(k*dt-.57)<10^(-8),
        vec=[v(k+1)-v(k-10); n(k+1)-n(k-10)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-b')
    end;
    if abs(k*dt-2.5)<10^(-8),
        vec=[v(k+1)-v(k-10); n(k+1)-n(k-10)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-b')
    end;
    if abs(k*dt-6.5)<10^(-8),
        vec=[v(k+1)-v(k-10); n(k+1)-n(k-10)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-b')
    end;
    if abs(k*dt-12.75)<10^(-8),
        vec=[v(k+1)-v(k-10); n(k+1)-n(k-10)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-b')
    end;
    if abs(k*dt-13.2)<10^(-8),
        vec=[v(k+1)-v(k-10); n(k+1)-n(k-10)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-b')
    end;
end;

plot(v,n,'-b','Linewidth',2);
hold off;

hold off;

shg;


