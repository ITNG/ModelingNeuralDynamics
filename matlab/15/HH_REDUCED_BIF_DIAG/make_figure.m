clear; clf;
tic
c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;

f=@(v) g_na*m_inf(v).^3.*(0.83-n_inf(v)).*(v_na-v)+ ...
       g_k*n_inf(v).^4.*(v_k-v)+ ...
       g_l*(v_l-v);

subplot(111);
t_final=300;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
v=zeros(m_steps+1,1);
n=zeros(m_steps+1,1);
m=zeros(m_steps+1,1);
h=zeros(m_steps+1,1);
n_sc=0;
n_uc=0;

alpha=5; beta=10;
i_ext_vec=alpha+[0:100]/100*(beta-alpha);
green_i=[]; green_v=[];
red_i=[]; red_v=[];
for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    v_left=-100;
    v_right=50;
    while v_right-v_left>10^(-10),
        v_c=(v_left+v_right)/2;
        if (f(v_c)+i_ext)*(f(v_left)+i_ext)>0,
            v_left=v_c;
        else
            v_right=v_c;
        end;
    end;
    v_c=(v_left+v_right)/2;
    fp_vec(ijk)=v_c;
    n_c=n_inf(v_c);
    J(1,1,ijk)=g_na*3*m_inf(v_c)^2*m_inf_p(v_c)*...
        (0.83-n_c)*(v_na-v_c)- ...
        g_na*m_inf(v_c)^3*(0.83-n_c)- ...
        g_k*n_c^4-g_l;
    J(1,2,ijk)=-g_na*m_inf(v_c)^3*(v_na-v_c)+...
        4*g_k*n_c^3*(v_k-v_c);
    J(2,1,ijk)=alpha_n_p(v_c)*(1-n_c)-beta_n_p(v_c)*n_c;
    J(2,2,ijk)=-alpha_n(v_c)-beta_n(v_c);
    E=eig(J(:,:,ijk));
    if real(E(1))<0 & real(E(2))<0 & abs(imag(E(1)))>10^(-4),
        red_v=[red_v,v_c];
        red_i=[red_i,i_ext];
        v(1)=v_c*1.01;      % start near the stable fixed point,
        n(1)=n_c*0.99;      % but not quite there, to find the
                            % unstable cycle
        m(1)=m_inf(v(1));
        h(1)=0.83-n(1);
        for k=1:m_steps,
            v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
                g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
            n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
            v_tmp=v(k)-dt05*v_inc;
            n_tmp=n(k)-dt05*n_inc;
            m_tmp=m_inf(v_tmp);
            h_tmp=0.83-n_tmp;
            v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
                g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
            n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
            v(k+1)=v(k)-dt*v_inc;
            n(k+1)=n(k)-dt*n_inc;
            m(k+1)=m_inf(v(k+1));
            h(k+1)=0.83-n(k+1);
        end;
        Bmax=max(v(round(2/3*m_steps:m_steps+1)));
        Bmin=min(v(round(2/3*m_steps:m_steps+1)));
        if Bmax<200,
            n_uc=n_uc+1;
            i_ext_uc(n_uc)=i_ext;
            Bmax_vec(n_uc)=Bmax;
            Bmin_vec(n_uc)=Bmin;
            v_c_uc(n_uc)=v_c;
        end;
    end;

    if real(E(1))>0 & real(E(2))>0 & abs(imag(E(1)))>10^(-4),
        green_v=[green_v,v_c];
        green_i=[green_i,i_ext];
    end;
    imag_part(ijk)=imag(E(1));
    real_part(ijk)=real(E(1));
    v_c_vec(ijk)=v_c;
    
    
    v(1)=70;       % start far away to get to the stable cycle
    n(1)=0.5;      
                     
    m(1)=m_inf(v(1));
    h(1)=0.83-n(1);
    for k=1:m_steps,
        v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
            g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
        n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
        v_tmp=v(k)+dt05*v_inc;
        n_tmp=n(k)+dt05*n_inc;
        m_tmp=m_inf(v_tmp);
        h_tmp=0.83-n_tmp;
        v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
            g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
        n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
        v(k+1)=v(k)+dt*v_inc;
        n(k+1)=n(k)+dt*n_inc;
        m(k+1)=m_inf(v(k+1));
        h(k+1)=0.83-n(k+1);
    end;
    Amax=max(v(round(2/3*m_steps:m_steps+1)));
    Amin=min(v(round(2/3*m_steps:m_steps+1)));
    if Amax-Amin>2,
        n_sc=n_sc+1;
        i_ext_sc(n_sc)=i_ext;
        Amax_vec(n_sc)=Amax;
        Amin_vec(n_sc)=Amin;
    end;
end;
plot(green_i,green_v,'--g','Linewidth',4)
hold on;
plot(red_i,red_v,'-r','Linewidth',4)
plot(i_ext_sc,Amax_vec,'-k','Linewidth',1);
plot(i_ext_sc,Amin_vec,'-k','Linewidth',1);
i_ext_uc=[i_ext_sc(1),i_ext_uc];
Bmax_vec=[Amax_vec(1),Bmax_vec];
Bmin_vec=[Amin_vec(1),Bmin_vec];
plot(i_ext_uc,Bmax_vec,'--k','Linewidth',2);
plot(i_ext_uc,Bmin_vec,'--k','Linewidth',2);



set(gca,'Fontsize',24);
set(gca,'Xtick',[]);
ylabel('$v$','Fontsize',32);
axis([5,10,-90,50])
text(4.5,-100,'$I_\ast \approx 5.25$','Fontsize',24);
text(7.4-0.75,-100,'$I_c \approx 7.4$','Fontsize',24);
plot([5.25,5.25],[-92,-87],'-k','Linewidth',1);
plot([7.4,7.4],[-92,-87],'-k','Linewidth',1);
hold off;
axis('square');
shg;
toc
