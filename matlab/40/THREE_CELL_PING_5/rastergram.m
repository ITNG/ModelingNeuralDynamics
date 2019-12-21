subplot(311);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b','Markersize',30); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r','Markersize',30);  hold on; end; 
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',14); 
set(gca,'Ytick',[1,2,3]);

axis([0,t_final,0,num_e+num_i+1]); 

subplot(312);
t=[1:m_steps]*dt;
plot(t,g_12,'-k','Linewidth',2);
set(gca,'Fontsize',14);
xlabel('$t$ [ms]','Fontsize',18);

shg;

