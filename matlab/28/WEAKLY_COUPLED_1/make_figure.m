clf; clear;

phi_B_0=0.4;                    % A is started at phase 0, and B at
                                % this phase. 

epsilon=0.5;                    % This is the small factor in the definition
                                % of the phase response function.

dt=0.01; dt05=dt/2;             % Time step used to solve the approximate
                                % differential equation for the phase
                                % difference. 


phi_A=0; phi_B=phi_B_0;         % initial phases of the two neurons 
num_spikes_A=1; num_spikes_B=0; % num_spikes_X is the number of spikes
                                % of neuron X that have occurred --- X=A
                                % or X=B. The initial assumption is that A
                                % just fired, and has had its effect on B. 
t_spikes_A=[0]; t_spikes_B=[];  % These are the spike times. 

t_vec(1)=0;                     % t_vec contains the spike times of either 
                                % A or B. 
ij=1;                           % ij is the length of t_vec.

psi_vec=[];
psi_vec(1)=phi_B-phi_A;         % psi_vec(ij) is the phase difference right 
                                % after the spike at time t_vec(ij).
                                
t_final=6/epsilon;              % We simulate up to time t_final.
                              
t=0;

while t<t_final                    % Simulate the two oscillators, one  
   delta(1)=ceiling(phi_A)-phi_A;  % spike at a time.
   delta(2)=ceiling(phi_B)-phi_B;  % "ceiling" is a self-defined function
   if delta(2)<delta(1),           % which differs from ceil in that
       t=t+delta(2);               % ceiling(x)=x+1 if x is an integer.
       num_spikes_B=num_spikes_B+1;
       t_spikes_B(num_spikes_B)=t;
       phi_A=phi_A+delta(2);
       phi_A=phi_A+epsilon*g(phi_A);
       phi_B=ceiling(phi_B);
       ij=ij+1;
       t_vec(ij)=t;
       psi_vec(ij)=phi_B-phi_A;
   else
       t=t+delta(1);               
       num_spikes_A=num_spikes_A+1;
       t_spikes_A(num_spikes_A)=t;
       phi_B=phi_B+delta(1);
       phi_B=phi_B+epsilon*g(phi_B);
       phi_A=ceiling(phi_A);
       ij=ij+1;
       t_vec(ij)=t;
       psi_vec(ij)=phi_B-phi_A;
   end;
end;


subplot(211);
for n=1:ij-1,
    plot([t_vec(n),t_vec(n+1)],[psi_vec(n),psi_vec(n)],'-k', ...
        'Linewidth',2);             % Note: between spikes, the phase
    hold on;                        % difference remains constant.
end;


m_steps=round(t_final/dt);
psi=zeros(m_steps+1,1);
psi(1)=phi_B_0;
for k=1:m_steps,                    % Now solve the differential equation.
    psi_inc=epsilon* ...
        (g(psi(k))-g(-psi(k)));     
    psi_tmp=psi(k)+dt05*psi_inc;    
    psi_inc=epsilon*...
        (g(psi_tmp)-g(-psi_tmp));
    psi(k+1)=psi(k)+dt*psi_inc;
end;

hold on;
plot([0:m_steps]*dt,psi,'-r','Linewidth',2);
set(gca,'Fontsize',16);
ylabel('$\psi$','Fontsize',20);
epsilon_str=num2str(epsilon);
title(['$\epsilon=$',epsilon_str],'Fontsize',20);
m1=max(psi_vec); m2=max(psi); M=max(m1,m2); M=ceil(M);
m1=min(psi_vec); m2=min(psi); m=min(m1,m2); m=floor(m);
axis([0,t_final,m,M]);
hold off;


epsilon=0.1;                    % This is the small factor in the definition
                                % of the phase response function.


phi_A=0; phi_B=phi_B_0;         % initial phases of the two neurons 
num_spikes_A=1; num_spikes_B=0; % num_spikes_X is the number of spikes
                                % of neuron X that have occurred --- X=A
                                % or X=B. The initial assumption is that A
                                % just fired, and has had its effect on B. 
t_spikes_A=[0]; t_spikes_B=[];  % These are the spike times. 

t_vec(1)=0;                     % t_vec contains the spike times of either 
                                % A or B. 
ij=1;                           % ij is the length of t_vec.

psi_vec=[];
psi_vec(1)=phi_B-phi_A;         % psi_vec(ij) is the phase difference right 
                                % after the spike at time t_vec(ij).
                                
t_final=6/epsilon;              % We simulate up to time t_final.
                              
t=0;

while t<t_final                    % Simulate the two oscillators, one  
   delta(1)=ceiling(phi_A)-phi_A;  % spike at a time.
   delta(2)=ceiling(phi_B)-phi_B;  % "ceiling" is a self-defined function
   if delta(2)<delta(1),           % which differs from ceil in that
       t=t+delta(2);               % ceiling(x)=x+1 if x is an integer.
       num_spikes_B=num_spikes_B+1;
       t_spikes_B(num_spikes_B)=t;
       phi_A=phi_A+delta(2);
       phi_A=phi_A+epsilon*g(phi_A);
       phi_B=ceiling(phi_B);
       ij=ij+1;
       t_vec(ij)=t;
       psi_vec(ij)=phi_B-phi_A;
   else
       t=t+delta(1);               
       num_spikes_A=num_spikes_A+1;
       t_spikes_A(num_spikes_A)=t;
       phi_B=phi_B+delta(1);
       phi_B=phi_B+epsilon*g(phi_B);
       phi_A=ceiling(phi_A);
       ij=ij+1;
       t_vec(ij)=t;
       psi_vec(ij)=phi_B-phi_A;
   end;
end;


subplot(212);
for n=1:ij-1,
    plot([t_vec(n),t_vec(n+1)],[psi_vec(n),psi_vec(n)],'-k', ...
        'Linewidth',2);             % Note: between spikes, the phase
    hold on;                        % difference remains constant.
end;


m_steps=round(t_final/dt);
psi=zeros(m_steps+1,1);
psi(1)=phi_B_0;
for k=1:m_steps,                    % Now solve the differential equation.
    psi_inc=epsilon* ...
        (g(psi(k))-g(-psi(k)));     
    psi_tmp=psi(k)+dt05*psi_inc;    
    psi_inc=epsilon*...
        (g(psi_tmp)-g(-psi_tmp));
    psi(k+1)=psi(k)+dt*psi_inc;
end;

hold on;
plot([0:m_steps]*dt,psi,'-r','Linewidth',2);
set(gca,'Fontsize',16);
ylabel('$\psi$','Fontsize',20);
xlabel('$t$ [units of $T$]','Fontsize',20);
epsilon_str=num2str(epsilon);
title(['$\epsilon=$',epsilon_str],'Fontsize',20);
axis([0,t_final,0,1]);
hold off;
m1=max(psi_vec); m2=max(psi); M=max(m1,m2); M=ceil(M);
m1=min(psi_vec); m2=min(psi); m=min(m1,m2); m=floor(m);
axis([0,t_final,m,M]);
shg;
