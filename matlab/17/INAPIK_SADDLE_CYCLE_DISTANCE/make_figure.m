clear; clf;

c=1;
g_na=20;
g_k=10; 
g_l=8;
v_na=60;
v_k=-90;
v_l=-80;
tau_n=0.15;

low=-1.38;
high=0;
i_ext_vec=low+[0:100]/100*(high-low);

for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);

    dt=0.001; dt05=dt/2;

    t_final=20;
    m_steps=round(t_final/dt);
    v=zeros(m_steps+1,1);
    m=zeros(m_steps+1,1);
    n=zeros(m_steps+1,1);

    v(1)=-30; 
    m(1)=m_inf(v(1));
    n(1)=0.2; 

    for k=1:m_steps,

        v_inc=(g_na*m(k)*(v_na-v(k))+ ...
            g_k*n(k)*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
        n_inc=(n_inf(v(k))-n(k))/tau_n;

        v_tmp=v(k)+dt05*v_inc;
        m_tmp=m_inf(v_tmp);
        n_tmp=n(k)+dt05*n_inc;

        v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
            g_k*n_tmp*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
        n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;

        v(k+1)=v(k)+dt*v_inc;
        m(k+1)=m_inf(v(k+1));
        n(k+1)=n(k)+dt*n_inc;
    end;


    f=@(v) g_na*m_inf(v).*(v_na-v)+g_k*n_inf(v).*(v_k-v)+g_l*(v_l-v);

    w=-100+[0:10000]/10000*150;

    k=[1:length(w)-1];
    ind=find((f(w(k))+i_ext).*(f(w(k+1))+i_ext)<=0);
    for klm=1:length(ind),
        j=ind(klm);
        w_low=w(j); w_high=w(j+1);
        while w_high-w_low>10^(-12),
            w_c=(w_low+w_high)/2;
            if (f(w_c)+i_ext)*(f(w_high)+i_ext)<=0,
                w_low=w_c;
            else
                w_high=w_c;
            end;
        end;
        v_c=(w_low+w_high)/2;
        n_c=n_inf(v_c);
        J=zeros(2,2);
        J(1,1)=g_na*m_inf_p(v_c)*(v_na-v_c)-g_na*m_inf(v_c)-g_k*n_c-g_l;
        J(1,2)=g_k*(v_k-v_c);
        J(2,1)=n_inf_p(v_c)/tau_n;
        J(2,2)=-1/tau_n;
        E=eig(J);
        if abs(imag(E(1)))<10^(-12) & real(E(1))*real(E(2))<0,
            v_star=v_c;
            n_star=n_c;
        end;

    end;

    d=min((v-v_star).^2+(n-n_star).^2);
    d_vec(ijk)=sqrt(d);
end;

subplot(211);
plot(i_ext_vec,d_vec,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$d$','Fontsize',20);
shg;

