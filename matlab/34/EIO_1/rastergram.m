subplot(311);

if num_spikes_o>0, plot(t_o_spikes,i_o_spikes,'.g');  hold on; end; 
if num_spikes_i>0, plot(t_i_spikes,i_i_spikes+num_o,'.b'); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i+num_o,'.r');  hold on; end; 

plot([0,t_final],[num_o+1/2,num_o+1/2],'--k','Linewidth',1);
plot([0,t_final],[num_i+num_o+1/2,num_i+num_o+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[num_o,num_i+num_o,num_e+num_i+num_o]);

axis([0,t_final,0,num_e+num_i+num_o+1]); 

subplot(312)

plot((0:m_steps)*dt,lfp_v,'-k','Linewidth',2);
set(gca,'Fontsize',16); 
axis([0,t_final,-90,-40]);


ylabel('mean($v_E$)','Fontsize',20);

subplot(313)

plot((0:m_steps)*dt,lfp_s,'-k','Linewidth',2);
set(gca,'Fontsize',16); 
axis([0,t_final,0,0.3]);

xlabel('$t$ [ms]','Fontsize',20);
ylabel('mean($s_E$)','Fontsize',20);
shg;

