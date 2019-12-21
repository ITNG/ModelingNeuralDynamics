subplot(211);

if num_spikes_e>0, 
    plot(t_e_spikes,i_e_spikes+num_i,'.r','Markersize',20);   
end; 
set(gca,'Fontsize',16); 
set(gca,'Ytick',[51,70]);
axis([0,t_final,50.5,70.5]); 
xlabel('$t$ [ms]','Fontsize',20);
hold on;

t=(0:m_steps)'*dt;
rate=zeros(m_steps+1,1);
sigma=3;
for k=1:num_spikes_e,
    rate=rate+exp(-(t-t_e_spikes(k)).^2/(2*sigma^2));
end;

% "rate" cannot be thought of as the instantaneous firing rate, but 
% as propotional to it, and since here we only want to find local maxima 
% of the instantaneous firing rate, that is good enough. 

rate_c=rate(2:m_steps);
rate_l=rate(1:m_steps-1);
rate_r=rate(3:m_steps+1);
ind=find(rate_c>rate_l & rate_c>rate_r);
t_vec=(ind-1)*dt;
for k=1:length(t_vec),
    plot([t_vec(k),t_vec(k)],[50.5,70.5],'-k','Linewidth',1);
end;
hold off;


shg;

