subplot(211);

plot(t_e_spikes,i_e_spikes,'.r','Markersize',20);  
set(gca,'Fontsize',16); 
xlabel('$t$ [ms]','Fontsize',20);
axis([50,t_final,71.5, 78.5]);
ylabel('E-cell number','Fontsize',20);




shg;

