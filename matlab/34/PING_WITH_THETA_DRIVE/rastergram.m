subplot(311);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b'); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r');  hold on; end; 
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',14); 
set(gca,'Ytick',[num_i,num_e+num_i]);

axis([0,t_final,0,num_e+num_i+1]); 

subplot(312);
plot((0:m_steps)*dt,lfp_v,'-k','Linewidth',2);
set(gca,'Fontsize',14);
title('mean($v_E$)','Fontsize',18);
axis([0,t_final,min(lfp_v)-10,max(lfp_v)+10]);

subplot(313);
plot((0:m_steps)*dt,lfp_s,'-k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$t$ [ms]','Fontsize',18);
title('mean($s_E$)','Fontsize',18);
axis([0,t_final,0,0.3]);

shg;

