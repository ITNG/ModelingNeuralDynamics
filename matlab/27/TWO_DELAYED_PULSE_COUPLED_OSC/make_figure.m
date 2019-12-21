clf; clear;

delta=0.1;
phi_A=0; phi_B=0.9;             % initial phases of the two neurons 
t_A_to_B=delta;                 % time in which next input reaches B
                                % Time is measured in units of T here.
t_B_to_A=inf;                   % time in which next input reaches A
t_present=0;                    % present time
t_final=20;                     % time to which to compute

num_spikes_A=1; num_spikes_B=0; % num_spikes_X is the number of spikes
                                % of neuron X that have occurred --- X=A
                                % or X=B. The initial assumption is that A
                                % just fired, and has had its effect on B. 
t_spikes_A=[0]; t_spikes_B=[];  % These are the spike times. 

                                

while t_present<t_final,
   T_vec=[1-phi_A,1-phi_B,t_A_to_B,t_B_to_A]; 
   T_0=min(T_vec);              % Time to next event of interest.
   done=0;
   
   if T_0==1-phi_A,
       phi_B=phi_B+1-phi_A;
       t_B_to_A=t_B_to_A-(1-phi_A);
       t_A_to_B=delta;
       t_present=t_present+1-phi_A;
       num_spikes_A=num_spikes_A+1;
       t_spikes_A(num_spikes_A)=t_present;
       phi_A=0;
       done=1;
   end;
   
   if T_0==1-phi_B & done==0,
       phi_A=phi_A+1-phi_B;
       t_A_to_B=t_A_to_B-(1-phi_B);
       t_B_to_A=delta;
       t_present=t_present+1-phi_B;
       num_spikes_B=num_spikes_B+1;
       t_spikes_B(num_spikes_B)=t_present;
       phi_B=0;
       done=1;
   end;
   
   if T_0==t_A_to_B & done==0,
       phi_B=f(phi_B+t_A_to_B);
       phi_A=phi_A+t_A_to_B;
       t_B_to_A=t_B_to_A-t_A_to_B;
       t_present=t_present+t_A_to_B;
       t_A_to_B=inf;
       done=1;
   end;
   
   if T_0==t_B_to_A & done==0,
       phi_A=f(phi_A+t_B_to_A);
       phi_B=phi_B+t_B_to_A;
       t_A_to_B=t_A_to_B-t_B_to_A;
       t_present=t_present+t_B_to_A;
       t_B_to_A=inf;
       done=1;
   end;
       
end;
t_spikes_A(num_spikes_A-2)
v=[
t_spikes_B(num_spikes_B-2)
t_spikes_B(num_spikes_B-1)
t_spikes_B(num_spikes_B-3)];
d=abs(v-t_spikes_A(num_spikes_A-2));
ind=find(d==min(d));
v(ind)

subplot(211);
plot(t_spikes_A,ones(num_spikes_A,1),'.r','Markersize',20);
hold on;
plot(t_spikes_B,2*ones(num_spikes_B,1),'.g','Markersize',20);
axis([t_final-10,t_final,0,3]);
set(gca,'Ytick',[]);
set(gca,'Fontsize',16);
%xlabel('$t$ [units of $T$]','Fontsize',20);
delta_str=num2str(delta);
title(['$\delta=$',delta_str,':'],'Fontsize',20);
hold off;

delta=0.7;
phi_A=0; phi_B=0.9;             % initial phases of the two neurons 
t_A_to_B=delta;                 % time in which next input reaches B
                                % Time is measured in units of T here.
t_B_to_A=inf;                   % time in which next input reaches A
t_present=0;                    % present time

num_spikes_A=1; num_spikes_B=0; % num_spikes_X is the number of spikes
                                % of neuron X that have occurred --- X=A
                                % or X=B. The initial assumption is that A
                                % just fired, and has had its effect on B. 
t_spikes_A=[0]; t_spikes_B=[];  % These are the spike times. 

                                

while t_present<t_final,
   T_vec=[1-phi_A,1-phi_B,t_A_to_B,t_B_to_A]; 
   T_0=min(T_vec);              % Time to next event of interest.
   done=0;
   
   if T_0==1-phi_A,
       phi_B=phi_B+1-phi_A;
       t_B_to_A=t_B_to_A-(1-phi_A);
       t_A_to_B=delta;
       t_present=t_present+1-phi_A;
       num_spikes_A=num_spikes_A+1;
       t_spikes_A(num_spikes_A)=t_present;
       phi_A=0;
       done=1;
   end;
   
   if T_0==1-phi_B & done==0,
       phi_A=phi_A+1-phi_B;
       t_A_to_B=t_A_to_B-(1-phi_B);
       t_B_to_A=delta;
       t_present=t_present+1-phi_B;
       num_spikes_B=num_spikes_B+1;
       t_spikes_B(num_spikes_B)=t_present;
       phi_B=0;
       done=1;
   end;
   
   if T_0==t_A_to_B & done==0,
       phi_B=f(phi_B+t_A_to_B);
       phi_A=phi_A+t_A_to_B;
       t_B_to_A=t_B_to_A-t_A_to_B;
       t_present=t_present+t_A_to_B;
       t_A_to_B=inf;
       done=1;
   end;
   
   if T_0==t_B_to_A & done==0,
       phi_A=f(phi_A+t_B_to_A);
       phi_B=phi_B+t_B_to_A;
       t_A_to_B=t_A_to_B-t_B_to_A;
       t_present=t_present+t_B_to_A;
       t_B_to_A=inf;
       done=1;
   end;
       
end;
t_spikes_A(num_spikes_A-2)
v=[
t_spikes_B(num_spikes_B-2)
t_spikes_B(num_spikes_B-1)
t_spikes_B(num_spikes_B-3)];
d=abs(v-t_spikes_A(num_spikes_A-2));
ind=find(d==min(d));
v(ind)

subplot(212);
plot(t_spikes_A,ones(num_spikes_A,1),'.r','Markersize',20);
hold on;
plot(t_spikes_B,2*ones(num_spikes_B,1),'.g','Markersize',20);
axis([t_final-10,t_final,0,3]);
set(gca,'Ytick',[]);
set(gca,'Fontsize',16);
xlabel('$t$ [units of $T$]','Fontsize',20);
delta_str=num2str(delta);
title(['$\delta=$',delta_str,':'],'Fontsize',20);
hold off;
shg;