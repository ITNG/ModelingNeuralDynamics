subplot(211);
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes,'.r');  hold on; end; 
hold off;
set(gca,'Fontsize',16); 
axis([0,t_final,0,num_e+1]); 
xlabel('$t$ [ms]','Fontsize',20);


shg;

