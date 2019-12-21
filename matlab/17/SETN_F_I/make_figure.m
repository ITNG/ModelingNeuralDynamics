clear; clf;

% This is very similar to the code in HH_F_I_CURVE. See there for 
% more extensive comments. 

tau_z=100;
z_max=0.05;

I_L=-0.06; I_R=0.02;
i_ext_vec=I_L+[0:30]/30*(I_R-I_L);

dt=0.01; dt05=dt/2;
max_num_spikes=3; t_max=1000; 
t_spikes=zeros(max_num_spikes,1);
z=zeros(round(t_max/dt)+1,1);
theta=z;

for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    theta(1)=0; 
    z(1)=0;
    k=1;
    t=0;
    num_spikes=0;
    while num_spikes<max_num_spikes && t < t_max,
        theta_inc=1-cos(theta(k))+(i_ext+z(k))*(1+cos(theta(k)));
        z_inc=-z(k)/tau_z+10*exp(-5*(1+cos(theta(k))))*(z_max-z(k));
        theta_tmp=theta(k)+dt05*theta_inc;
        z_tmp=z(k)+dt05*z_inc;
        theta_inc=1-cos(theta_tmp)+(i_ext+z_tmp)*(1+cos(theta_tmp));
        z_inc=-z_tmp/tau_z+10*exp(-5*(1+cos(theta_tmp)))*(z_max-z_tmp);
        theta(k+1)=theta(k)+dt*theta_inc;
        z(k+1)=z(k)+dt*z_inc;
        if theta(k+1)>pi,
            num_spikes=num_spikes+1;
            t_spikes(num_spikes)= ...
                (k-1)*dt*(theta(k+1)-pi)+k*dt*(pi-theta(k));
            t_spikes(num_spikes)= ...
                 t_spikes(num_spikes)/(theta(k+1)-theta(k));
            theta(k+1)=theta(k+1)-2*pi;
        end;
        k=k+1; t=t+dt;
    end;
    if num_spikes==max_num_spikes,
        f_low(ijk)=1000/ ...
            (t_spikes(max_num_spikes)-t_spikes(max_num_spikes-1));
    else
        f_low(ijk)=0;
    end;
end;

for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    theta(1)=9/10*pi; 
    z(1)=0;
    k=1;
    t=0;
    num_spikes=0;
    while num_spikes<max_num_spikes && t < t_max,
        theta_inc=1-cos(theta(k))+(i_ext+z(k))*(1+cos(theta(k)));
        z_inc=-z(k)/tau_z+10*exp(-5*(1+cos(theta(k))))*(z_max-z(k));
        theta_tmp=theta(k)+dt05*theta_inc;
        z_tmp=z(k)+dt05*z_inc;
        theta_inc=1-cos(theta_tmp)+(i_ext+z_tmp)*(1+cos(theta_tmp));
        z_inc=-z_tmp/tau_z+10*exp(-5*(1+cos(theta_tmp)))*(z_max-z_tmp);
        theta(k+1)=theta(k)+dt*theta_inc;
        z(k+1)=z(k)+dt*z_inc;
        if theta(k+1)>pi,
            num_spikes=num_spikes+1;
            t_spikes(num_spikes)= ...
                (k-1)*dt*(theta(k+1)-pi)+k*dt*(pi-theta(k));
            t_spikes(num_spikes)= ...
                 t_spikes(num_spikes)/(theta(k+1)-theta(k));
            theta(k+1)=theta(k+1)-2*pi;
        end;
        k=k+1; t=t+dt;
    end;
    if num_spikes==max_num_spikes,
        f_high(ijk)=1000/ ...
            (t_spikes(max_num_spikes)-t_spikes(max_num_spikes-1));
    else
        f_high(ijk)=0;
    end;
end;

subplot(211);
plot(i_ext_vec,f_low,'.k','Markersize',15);
hold on;
plot(i_ext_vec,f_high,'ok','Markersize',10,'Linewidth',1);
axis([min(i_ext_vec),max(i_ext_vec),0,100]);
set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$f$','Fontsize',20);
set(gca,'Xtick',[-0.06:0.02:0.02]);
shg;
hold off;
