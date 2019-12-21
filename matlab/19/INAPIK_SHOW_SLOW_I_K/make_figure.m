clear; clf;

c=1;
g_na=20;
g_k=10; 
g_l=8;
v_na=60;
v_k=-90;
v_l=-80;
tau_n=0.15;

i_ext=7;

g_k_slow=5; tau_n_slow=20;


t_final=100; dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); v=z; m=z; n=z; n_slow=z;
v(1)=-70; m(1)=m_inf(v(1)); n(1)=0.6; n_slow(1)=0;

for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+ ...
        g_k_slow*n_slow(k)*(v_k-v(k))+ ...
        g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    n_slow_inc=(n_slow_inf(v(k))-n_slow(k))/tau_n_slow;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    n_slow_tmp=n_slow(k)+dt05*n_slow_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+ ...
        g_k_slow*n_slow_tmp*(v_k-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    n_slow_inc=(n_slow_inf(v_tmp)-n_slow_tmp)/tau_n_slow;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
    n_slow(k+1)=n_slow(k)+dt*n_slow_inc;
    
end;

subplot(211);
t=[0:m_steps]*dt;
plot(t,v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
ylabel('$v$ [mV]','Fontsize',20);
axis([40,90,-80,0]);   % Which window to show, here and in the other
                       % subplot, was of course decided simply based on
                       % looking at the plots and finding a window that 
                       % would give a complete look at one burst and the
                       % interburst interval following it.

subplot(212);
plot(t,i_ext+g_k_slow*n_slow.*(v_k-v),'-k','Linewidth',2);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$I_{\rm eff}$ [$\mu$A/cm$^2$]','Fontsize',20);
hold on;
plot([40,90],[4.5,4.5],'-r','Linewidth',2);     % 4.5 and -1.4 are roughly
plot([40,90],[-1.4,-1.4],'-b','Linewidth',2);   % the values of I_c and I_star,
                                                % respectively, that is,
                                                % the critical values for
                                                % the two bifurcations.
                                                % Great precision would not
                                                % be useful here.
                                                
hold off;
set(gca,'Fontsize',16);
axis([40,90,-2,7]);
shg;