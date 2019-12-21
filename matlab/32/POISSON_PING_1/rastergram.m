subplot(211);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b'); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r');  hold on; end; 
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[num_i,num_e+num_i]);

axis([0,t_final,0,num_e+num_i+1]); 

subplot(212);
for k=round(5/dt):m_steps-round(5/dt),
    f_e(k-round(5/dt)+1)=sum(spike_count(k-round(5/dt)+1:k+round(5/dt)))/10;
end;
f_e=f_e/num_e*1000;
t_e=(round(5/dt):m_steps-round(5/dt))*dt-dt/2;
plot(t_e,f_e,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$f_E$ [Hz]','Fontsize',20);

shg;

