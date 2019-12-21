subplot(211);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b'); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r');  hold on; end; 
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[num_i,num_e+num_i]);

axis([0,t_final,0,num_e+num_i+1]); 

hold on;
t=[0:1000]/1000*t_final;
y=30*exp(cos(pi*(t/P-phi)).^8)+30;
plot(t,y,'-g','Linewidth',2);
hold off;



shg;

