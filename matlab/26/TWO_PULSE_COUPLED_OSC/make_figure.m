clf; clear;

phi_A=0; phi_B=0.5;             % initial phases of the two neurons 
num_spikes_A=1; num_spikes_B=0; % num_spikes_X is the number of spikes
                                % of neuron X that have occurred --- X=A
                                % or X=B. The initial assumption is that A
                                % just fired, and has had its effect on B. 
t_spikes_A=[0]; t_spikes_B=[];  % These are the spike times. 

N=12;                           % Compute until a total of N spikes have 
                                % occurred. 
t=0;

for k=1:N,
   t=t+(1-phi_B);               % The next event of interest is a spike
                                % of neuron B, which will occur in time
                                % 1 - phi_B. (The time unit is the
                                % intrinsic period of the two neurons.) 
   num_spikes_B=num_spikes_B+1;
   t_spikes_B(num_spikes_B)=t;
   phi_A=f(1-phi_B);            % Neuron A is at phase 1- phi_B just 
                                % before neuron B fires, and is re-set to
                                % phase f(1-phi_B). 
   phi_B=0;                     % Now neuron B is at phase 0. 
   t=t+(1-phi_A);               % All that just happened repeats itself, 
                                % with the roles of A and B reversed. 
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