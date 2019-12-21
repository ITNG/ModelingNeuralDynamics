subplot(212);
if num_spikes_i>0, plot(t_i_spikes,ones(num_spikes_i,1),'.b','Markersize',30); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,2*ones(num_spikes_e,1),'.r','Markersize',30);  end;
tt=0;
while tt<=t_final,
    plot([tt,tt],[0,3],'--k','Linewidth',1);
    tt=tt+T_main;
end;
hold off;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[]);
xlabel('$t$ [ms]','Fontsize',20);

axis([0,t_final,0,3]); 

subplot(211);
plot(t,I_main,'-k','Linewidth',2);
hold on;
plot(t,I_dist,'--k','Linewidth',2);
hold off;
set(gca,'Fontsize',16);



shg;

