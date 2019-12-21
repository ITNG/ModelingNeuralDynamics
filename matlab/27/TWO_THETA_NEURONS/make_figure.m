clear; clf;

subplot(111);
epsilon_vec=[1:1000]/1000*5;
delta_vec=1/2-1/pi*atan(epsilon_vec/2);
plot(epsilon_vec,delta_vec,'-k','Linewidth',5);
set(gca,'Fontsize',24);
xlabel('$\epsilon$','Fontsize',32);
ylabel('$\delta$','Fontsize',32);
axis([0,5,0,1]); axis('square');
shg;

h=patch([0,epsilon_vec,5,epsilon_vec(1000:-1:1),0,0], ...
        [0.5,delta_vec,1,ones(1,1000),1,0.5],'y','FaceAlpha',0.2);

    
tau_m=2;
I=0.14;
T=pi*tau_m/sqrt(tau_m*I-0.25);
t_final=10000;
dt=0.01;
m_steps=round(t_final/dt);
dt05=dt/2;

N=10;
for i=1:N-1,
    for j=1:N-1,
        delta=j/N;
        epsilon=i/N*5;
        dv=sqrt(tau_m*I-0.25)*epsilon;
        
        theta_A(1)=0;
        theta_B(1)=pi/6;
        num_spikes_A=0;
        num_spikes_B=0;
        t_spikes_A=[];
        t_spikes_B=[];

        for k=1:m_steps,

            ind=find(t_spikes_B+delta*T>(k-1)*dt&t_spikes_B+delta*T<=k*dt);
            if length(ind)>0,
                t_0=t_spikes_B(ind)+delta*T;
                dt_1=t_0-(k-1)*dt; dt_105=dt_1/2;
                theta_A_inc=-cos(theta_A(k))/tau_m+2*I*(1+cos(theta_A(k)));
                theta_A_tmp=theta_A(k)+dt_105*theta_A_inc;
                theta_A_inc=-cos(theta_A_tmp)/tau_m+2*I*(1+cos(theta_A_tmp));
                t_A=theta_A(k)+dt_1*theta_A_inc;
                t_A=2*atan(tan(t_A/2)+2*dv);
                dt_2=k*dt-t_0; dt_205=dt_2/2;
                theta_A_inc=-cos(t_A)/tau_m+2*I*(1+cos(t_A));
                theta_A_tmp=t_A+dt_205*theta_A_inc;
                theta_A_inc=-cos(theta_A_tmp)/tau_m+2*I*(1+cos(theta_A_tmp));
                theta_A(k+1)=t_A+dt_2*theta_A_inc;
            else
                theta_A_inc=-cos(theta_A(k))/tau_m+2*I*(1+cos(theta_A(k)));
                theta_A_tmp=theta_A(k)+dt05*theta_A_inc;
                theta_A_inc=-cos(theta_A_tmp)/tau_m+2*I*(1+cos(theta_A_tmp));
                theta_A(k+1)=theta_A(k)+dt*theta_A_inc;
            end;

            ind=find(t_spikes_A+delta*T>(k-1)*dt&t_spikes_A+delta*T<=k*dt);
            if length(ind)>0,
                t_0=t_spikes_A(ind)+delta*T;
                dt_1=t_0-(k-1)*dt; dt_105=dt_1/2;
                theta_B_inc=-cos(theta_B(k))/tau_m+2*I*(1+cos(theta_B(k)));
                theta_B_tmp=theta_B(k)+dt_105*theta_B_inc;
                theta_B_inc=-cos(theta_B_tmp)/tau_m+2*I*(1+cos(theta_B_tmp));
                t_B=theta_B(k)+dt_1*theta_B_inc;
                t_B=2*atan(tan(t_B/2)+2*dv);
                dt_2=k*dt-t_0; dt_205=dt_2/2;
                theta_B_inc=-cos(t_B)/tau_m+2*I*(1+cos(t_B));
                theta_B_tmp=t_B+dt_205*theta_B_inc;
                theta_B_inc=-cos(theta_B_tmp)/tau_m+2*I*(1+cos(theta_B_tmp));
                theta_B(k+1)=t_B+dt_2*theta_B_inc;
            else
                theta_B_inc=-cos(theta_B(k))/tau_m+2*I*(1+cos(theta_B(k)));
                theta_B_tmp=theta_B(k)+dt05*theta_B_inc;
                theta_B_inc=-cos(theta_B_tmp)/tau_m+2*I*(1+cos(theta_B_tmp));
                theta_B(k+1)=theta_B(k)+dt*theta_A_inc;
            end;

            if theta_A(k+1)>pi,
                num_spikes_A=num_spikes_A+1;
                tt=(k-1)*dt*(theta_A(k+1)-pi)+k*dt*(pi-theta_A(k));
                tt=tt/(theta_A(k+1)-theta_A(k));
                t_spikes_A(num_spikes_A)=tt;
                theta_A(k+1)=theta_A(k+1)-2*pi;
            end;

            if theta_B(k+1)>pi,
                num_spikes_B=num_spikes_B+1;
                tt=(k-1)*dt*(theta_B(k+1)-pi)+k*dt*(pi-theta_B(k));
                tt=tt/(theta_B(k+1)-theta_B(k));
                t_spikes_B(num_spikes_B)=tt;
                theta_B(k+1)=theta_B(k+1)-2*pi;
            end;


        end;
        d=zeros(num_spikes_A,1);
        for k=1:num_spikes_A,
            d(k)=min(abs(t_spikes_A(k)-t_spikes_B));
        end;
        sync_measure=d(num_spikes_A-1);
        hold on;
        if sync_measure<10^(-2),
            plot(epsilon,delta,'.r','Markersize',25);
        else
            plot(epsilon,delta,'.b','Markersize',25);
        end;
    end;
end;
hold off;
shg;
