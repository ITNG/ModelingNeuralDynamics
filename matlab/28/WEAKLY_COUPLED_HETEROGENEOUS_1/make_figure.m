clf; clear;
 
c=0.08;                         % T_B=(1+epsilon c)T_A.

epsilon=0.1;                    % This is the small factor in the definition
                                % of the phase response function.
                                
T_B=1+epsilon*c;                % We take T_A to be 1.
                                
phi_B_0=0.5;                    % A is started at phase 0, and B at
                                % this phase. 

dt=0.1; dt05=dt/2;              % Time step used to solve the approximate
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

t_final=500;


m_steps=round(t_final/dt);
psi=zeros(m_steps+1,1);
psi(1)=phi_B_0;
for k=1:m_steps,                    % Solve the differential equation. 
    psi_inc=epsilon*(g(psi(k))-g(-psi(k))) ...
        -c*epsilon;                 
    psi_tmp=psi(k)+dt05*psi_inc;    
    psi_inc=epsilon*(g(psi_tmp)-g(-psi_tmp)) ...
        -c*epsilon;
    psi(k+1)=psi(k)+dt*psi_inc;
end;
subplot(211);
plot([0:m_steps]*dt,psi,'-r','Linewidth',2);
hold on;


t=0;

while t<=t_final,
    delta(1)=(ceiling(phi_B)-phi_B)*T_B;   % time to next spike of B
    delta(2)=ceiling(phi_A)-phi_A;         % time to next spike of A
    if delta(1)<delta(2),                  % ceiling(x)=ceil(x) if x is
                                           % not an integer, and =x+1 if x
                                           % is an integer.
        t=t+delta(1);
        num_spikes_B=num_spikes_B+1;
        t_spikes_B(num_spikes_B)=t;
        phi_A=phi_A+delta(1);
        phi_A=phi_A+epsilon*g(phi_A);
        phi_B=ceiling(phi_B);
        ij=ij+1;
        t_vec(ij)=t;
        psi_vec(ij)=phi_B-phi_A;
    else
        t=t+delta(2);
        num_spikes_A=num_spikes_A+1;
        t_spikes_A(num_spikes_A)=t;
        phi_B=phi_B+delta(2)/T_B;
        phi_B=phi_B+epsilon*g(phi_B);
        phi_A=ceiling(phi_A);
        ij=ij+1;
        t_vec(ij)=t;
        psi_vec(ij)=phi_B-phi_A;
    end;
end;


for n=1:ij-1,
    plot([t_vec(n),t_vec(n+1)],[psi_vec(n),psi_vec(n)],'-k', ...
        'Linewidth',2);             % Note: between spikes, the phase
                                    % difference remains constant.
end;
set(gca,'Fontsize',16);
ylabel('$\psi$','Fontsize',20);



hold off;


epsilon_str=num2str(epsilon);
c_str=num2str(c);
title(['$\epsilon=$',epsilon_str, ',  $c=$',c_str],'Fontsize',20);
hold off;

m1=max(psi_vec); m2=max(psi); M=max(m1,m2); M=ceil(M);
m1=min(psi_vec); m2=min(psi); m=min(m1,m2); m=floor(m);
axis([0,t_final,m,M]);
shg;

p1=num2str(psi_vec(length(psi_vec)));
disp(['actual phase difference:  ',p1]); 
p2=num2str(psi(length(psi)));
disp(['phase difference predicted based on DE:  ',p2]); 


clear;
 
c=0.12;                         % T_B=(1+epsilon c)T_A.

epsilon=0.1;                    % This is the small factor in the definition
                                % of the phase response function.
                                
T_B=1+epsilon*c;                % We take T_A to be 1.
                                
phi_B_0=0.5;                    % A is started at phase 0, and B at
                                % this phase. 

dt=0.1; dt05=dt/2;              % Time step used to solve the approximate
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

t_final=500;

m_steps=round(t_final/dt);
psi=zeros(m_steps+1,1);
psi(1)=phi_B_0;
for k=1:m_steps,                    % Solve the differential equation. 
    psi_inc=epsilon*(g(psi(k))-g(-psi(k))) ...
        -c*epsilon;                 
    psi_tmp=psi(k)+dt05*psi_inc;    
    psi_inc=epsilon*(g(psi_tmp)-g(-psi_tmp)) ...
        -c*epsilon;
    psi(k+1)=psi(k)+dt*psi_inc;
end;
subplot(212);
plot([0:m_steps]*dt,psi,'-r','Linewidth',2);
hold on;


t=0;

while t<=t_final,
    delta(1)=(ceiling(phi_B)-phi_B)*T_B;   % time to next spike of B
    delta(2)=ceiling(phi_A)-phi_A;         % time to next spike of A
    if delta(1)<delta(2),                  % ceiling(x)=ceil(x) if x is
                                           % not an integer, and =x+1 if x
                                           % is an integer.
        t=t+delta(1);
        num_spikes_B=num_spikes_B+1;
        t_spikes_B(num_spikes_B)=t;
        phi_A=phi_A+delta(1);
        phi_A=phi_A+epsilon*g(phi_A);
        phi_B=ceiling(phi_B);
        ij=ij+1;
        t_vec(ij)=t;
        psi_vec(ij)=phi_B-phi_A;
    else
        t=t+delta(2);
        num_spikes_A=num_spikes_A+1;
        t_spikes_A(num_spikes_A)=t;
        phi_B=phi_B+delta(2)/T_B;
        phi_B=phi_B+epsilon*g(phi_B);
        phi_A=ceiling(phi_A);
        ij=ij+1;
        t_vec(ij)=t;
        psi_vec(ij)=phi_B-phi_A;
    end;
end;


for n=1:ij-1,
    plot([t_vec(n),t_vec(n+1)],[psi_vec(n),psi_vec(n)],'-k', ...
        'Linewidth',2);             % Note: between spikes, the phase
                                    % difference remains constant.
end;
set(gca,'Fontsize',16);
ylabel('$\psi$','Fontsize',20);
xlabel('$t$ [units of $T_A$]','Fontsize',20);


hold off;


epsilon_str=num2str(epsilon);
c_str=num2str(c);
title(['$\epsilon=$',epsilon_str, ',  $c=$',c_str],'Fontsize',20);
hold off;

m1=max(psi_vec); m2=max(psi); M=max(m1,m2); M=ceil(M);
m1=min(psi_vec); m2=min(psi); m=min(m1,m2); m=floor(m);
axis([0,t_final,m,M]);
shg;

p1=num2str(psi_vec(length(psi_vec)));
disp(['actual phase difference:  ',p1]); 
p2=num2str(psi(length(psi)));
disp(['phase difference predicted based on DE:  ',p2]); 

