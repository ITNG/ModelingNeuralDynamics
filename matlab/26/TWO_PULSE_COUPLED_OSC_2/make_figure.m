clf; clear;

% For explanations, see the code in TWO_PULSE_COUPLED_OSC.

phi_A=0;
phi_B=0.05; 
num_spikes_A=1;
t_spikes_A=[0];
num_spikes_B=0;
t_spikes_B=[];

N=12;
t=0;

for k=1:N,
   t=t+(1-phi_B);
   num_spikes_B=num_spikes_B+1;
   t_spikes_B(num_spikes_B)=t;
   phi_A=f(1-phi_B);
   phi_B=0;
   t=t+(1-phi_A);
   num_spikes_A=num_spikes_A+1;
   t_spikes_A(num_spikes_A)=t;
   phi_B=f(1-phi_A);
   phi_A=0;
end;

subplot(211);
plot(t_spikes_A,ones(num_spikes_A,1),'.r','Markersize',20);
hold on;
plot(t_spikes_B,2*ones(num_spikes_B,1),'.b','Markersize',20);
axis([0,N-2,0,3]);
set(gca,'Xtick',[1:N-1]);
set(gca,'Ytick',[]);
set(gca,'Fontsize',16);
xlabel('$t$ [units of $T$]','Fontsize',20);
title('spike times of A (red) and B (blue)','Fontsize',20);
hold off;
shg;