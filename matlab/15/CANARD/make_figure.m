clear; clf;

tic
%a=100/21; tau_n=100/sqrt(3); % parameters in the FN-equation
a=5; tau_n=60;               % parameters in the FN-equation
 % the canard explosion happens in this interval;



t_final=10000;              % we find the amplitude of the limit cycle
                            % by solving from time 0 to time t_final ...
dt=0.01; dt05=dt/2;         % ... with time step dt
m_steps=round(t_final/dt);
v=zeros(m_steps+1,1); 
n=zeros(m_steps+1,1);
subplot(111);
amp_vec=[1,2,3,3.5,3.78];
for ijk=1:5,
    amp_target=amp_vec(ijk);
    I_left=-4.28; I_right=-4.25;
    while I_right-I_left>10^(-15),
        I=(I_left+I_right)/2;
        v(1)=0; n(1)=0.1;
        for k=1:m_steps,
            v_inc=v(k)-v(k)^3/3-n(k)+I;
            n_inc=(a*v(k)-n(k))/tau_n;
            v_tmp=v(k)+dt05*v_inc;
            n_tmp=n(k)+dt05*n_inc;
            v_inc=v_tmp-v_tmp^3/3-n_tmp+I;
            n_inc=(a*v_tmp-n_tmp)/tau_n;
            v(k+1)=v(k)+dt*v_inc;
            n(k+1)=n(k)+dt*n_inc;
        end;
        vv=v(4*m_steps/5:m_steps+1); 
        amp=max(vv)-min(vv);
        if amp>amp_target,
            I_right=I;
        else
            I_left=I;
        end
    end;
    I
    amp

    vv=v(m_steps/2:m_steps+1);
    nn=n(m_steps/2:m_steps+1);
    subplot(221);
    if amp_target==3.5,
        plot(vv,nn,'-b','Linewidth',6);
    else
        plot(vv,nn,'-k','Linewidth',1);  
    end;
    set(gca,'Fontsize',14);
    xlabel('$v$','Fontsize',18);
    ylabel('$n$','Fontsize',18);
    A=-3; B=2; C=-6; D=-2;
    axis([A,B,C,D]); axis('square');
    hold on;
    
    if ijk==1,
        subplot(222);
        L=length(vv);
        L_200=round(200/dt);
        plot((0:L_200)*dt, ...
              vv(L-L_200:L),'-k','Linewidth',2);
        set(gca,'Fontsize',14);
        xlabel('$t$','Fontsize',18);
        ylabel('$v$','Fontsize',18);
        I_str=num2str(I,7);
        title(['$I=$',I_str],'Fontsize',14);
        axis([0,200,-3,2]);
    end;
    if ijk==2,
        subplot(223);
        L=length(vv);
        L_200=round(200/dt);
        plot((0:L_200)*dt, ...
              vv(L-L_200:L),'-k','Linewidth',2);
        set(gca,'Fontsize',14);
        xlabel('$t$','Fontsize',18);
        ylabel('$v$','Fontsize',18);
        I_str=num2str(I,7);
        title(['$I=$',I_str],'Fontsize',14);
        axis([0,200,-3,2]);
    end;
    if ijk==5,
        subplot(224);
        L=length(vv);
        L_200=round(200/dt);
        plot((0:L_200)*dt, ...
              vv(L-L_200:L),'-k','Linewidth',2);
        set(gca,'Fontsize',14);
        xlabel('$t$','Fontsize',18);
        ylabel('$v$','Fontsize',18);
        I_str=num2str(I,7);
        title(['$I=$',I_str],'Fontsize',14);
        axis([0,200,-3,2]);
    end;
end;
subplot(221);
hold on;
v=A+(0:100)/100*(B-A);
n=v-v.^3/3+I;
plot(v,n,'-g','Linewidth',2);
n=a*v;
plot(v,n,'-r','Linewidth',2);
hold off;
shg;


