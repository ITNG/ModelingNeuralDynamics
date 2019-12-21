subplot(211);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b','Markersize',15); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r','Markersize',15);  hold on; end; 
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[50,70]);


axis([t_final-200,t_final,49.5,70]);

shg;

