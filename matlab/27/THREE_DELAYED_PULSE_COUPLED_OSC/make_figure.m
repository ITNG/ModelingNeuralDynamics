clear; clf; rng('default'); rng(63806);

N=3;           % number of oscillators

delta=.45*ones(N,N);  % delta(i,j) is the conductance delay from oscillator
                      % i to oscillator j, in units of T. We assume
                      % delta(i,j)<1.

t_final=200;     % Time to which simulation is run. The time unit is the 
                 % intrinsic period of the oscillators. 

t_0=0;                  % initital time
phi=0.9+rand(N,1)*0.1;  % initial phases 
num_spikes=0;   % total number of spikes recorded thus far

M=0;            % total number of signals currently traveling

t=zeros(N^2,1); % The first M entries in t indicate how close the
                % traveling signals are to their targets. The array t is 
                % initialized because this makes Matlab run much more
                % efficiently. 
                
k_to=zeros(N^2,1);   % The first M entries in k_to indicate where the currently
                     % traveling signals are headed. 
k_from=zeros(N^2,1); % The first M entries in k_from indicate where the currently
                     % traveling signals originated.
                
tic

while t_0<t_final,
    next_spike=min(1-phi);       % time to next spike
    if M==0,
        next_signal=inf;
    else
        next_signal=min(t(1:M)); % time to next signal arrival
    end;
    next=min(next_spike,next_signal);   % time to next event (either spike
                                        % or signal arrival)
    if next==next_signal,
        j0=find(t(1:M)==next_signal);   % find the index of the signal 
                                        % that will arrive the soonest
                                        
        j0=max(j0);                     % Break ties by choosing, among
                                        % possibly several signals that 
                                        % will arrive soonest, the one
                                        % with the largest index j. (This
                                        % is an entirely arbitrary choice
                                        % that won't affect the
                                        % simulation.)
        phi=phi+next_signal;                  
        t_0=t_0+next_signal;
        t=t-next_signal;
        phi(k_to(j0))=phi(k_to(j0))+g(phi(k_to(j0)))/(N-1);  
                                        % the phase response function
                                        % is scaled by N-1 to keep the
                                        % total amount of input per
                                        % oscillator limited as N
                                        % increases.
                                                             
                                                    
        t(j0:M-1)=t(j0+1:M);            % this eliminates signal j0, which
        k_to(j0:M-1)=k_to(j0+1:M);      % just reached its target.
        k_from(j0:M-1)=k_from(j0+1:M);
        M=M-1;
        
    else
        
        i0=find(1-phi==next_spike);     % find the index of the neuron
                                        % that will spike the soonest,
        i0=max(i0);                     % breaking ties in an arbitrary
                                        % way
        phi=phi+next_spike;
        t_0=t_0+next_spike;
        t(1:M)=t(1:M)-next_spike;
        phi(i0)=0;
        
        num_spikes=num_spikes+1;        % record that a spike occurred
        t_spikes(num_spikes)=t_0;       % t_spikes(n) is the time of the
                                        % n-th spike
        i_spikes(num_spikes)=i0;        % at time t=t_spikes(n), the 
                                        % neuron that spikes is neuron
                                        % number i_spikes(n)
        ind=[M+1:M+N-1];                % N-1 new signals will now be
                                        % traveling. These are labeled
                                        % M+1 through M+N-1 (M is not
                                        % yet updated at this point). 
        k_from(ind)=i0;
        k_to(ind)=[1:i0-1,i0+1:N];
        for l=1:length(ind),
            t(ind(l))=delta(k_from(ind(l)),k_to(ind(l)));
        end;
        M=M+N-1;
    end;
end;
    
subplot(211);
for k=1:num_spikes, 
    plot(t_spikes(k),i_spikes(k),...    % plot spike rastergram
    '.k','Markersize',20); 
    hold on;
end;
set(gca,'Fontsize',16);
N_half=round(N/2);                      % beautify plot a bit
if N_half>1&N_half<N,
    set(gca,'Ytick',[1,N_half,N]);
else
    set(gca,'Ytick',[1:N]);
end;
    
axis([t_final-10,t_final,0,N+1]);        % show only last 5 time units
if max(max(delta))==min(min(delta)),
    delta_str=num2str(delta(1,2));
    title(['$\delta=$',delta_str,':'],'Fontsize',20);
end;
%xlabel('t','Fontsize',20);
hold off;
                                        
toc     % Matlab computes execution time for statements bracketed by
        % "tic" and "toc". This can be useful when you play with the
        % program.


delta=.55*ones(N,N);   % delta(i,j) is the conductance delay from oscillator
                      % i to oscillator j, in units of T. We assume
                      % delta(i,j)<1.


t_0=0;                  % initital time
phi=0.9+rand(N,1)*0.1;  % initial phases 
num_spikes=0;   % total number of spikes recorded thus far

M=0;            % total number of signals currently traveling

t=zeros(N^2,1); % The first M entries in t indicate how close the
                % traveling signals are to their targets. The array t is 
                % initialized because this makes Matlab run much more
                % efficiently. 
                
k_to=zeros(N^2,1);   % The first M entries in k_to indicate where the currently
                     % traveling signals are headed. 
k_from=zeros(N^2,1); % The first M entries in k_from indicate where the currently
                     % traveling signals originated.
                
tic

while t_0<t_final,
    next_spike=min(1-phi);       % time to next spike
    if M==0,
        next_signal=inf;
    else
        next_signal=min(t(1:M)); % time to next signal arrival
    end;
    next=min(next_spike,next_signal);   % time to next event (either spike
                                        % or signal arrival)
    if next==next_signal,
        j0=find(t(1:M)==next_signal);   % find the index of the signal 
                                        % that will arrive the soonest
                                        
        j0=max(j0);                     % Break ties by choosing, among
                                        % possibly several signals that 
                                        % will arrive soonest, the one
                                        % with the largest index j. (This
                                        % is an entirely arbitrary choice
                                        % that won't affect the
                                        % simulation.)
        phi=phi+next_signal;                  
        t_0=t_0+next_signal;
        t=t-next_signal;
        phi(k_to(j0))=phi(k_to(j0))+g(phi(k_to(j0)));
                                          
                                                    
        t(j0:M-1)=t(j0+1:M);            % this eliminates signal j0, which
        k_to(j0:M-1)=k_to(j0+1:M);      % just reached its target.
        k_from(j0:M-1)=k_from(j0+1:M);
        M=M-1;
        
    else
        
        i0=find(1-phi==next_spike);     % find the index of the neuron
                                        % that will spike the soonest,
        i0=max(i0);                     % breaking ties in an arbitrary
                                        % way
        phi=phi+next_spike;
        t_0=t_0+next_spike;
        t(1:M)=t(1:M)-next_spike;
        phi(i0)=0;
        
        num_spikes=num_spikes+1;        % record that a spike occurred
        t_spikes(num_spikes)=t_0;       % t_spikes(n) is the time of the
                                        % n-th spike
        i_spikes(num_spikes)=i0;        % at time t=t_spikes(n), the 
                                        % neuron that spikes is neuron
                                        % number i_spikes(n)
        ind=[M+1:M+N-1];                % N-1 new signals will now be
                                        % traveling. These are labeled
                                        % M+1 through M+N-1 (M is not
                                        % yet updated at this point). 
        k_from(ind)=i0;
        k_to(ind)=[1:i0-1,i0+1:N];
        for l=1:length(ind),
            t(ind(l))=delta(k_from(ind(l)),k_to(ind(l)));
        end;
        M=M+N-1;
    end;
end;
    
subplot(212);
for k=1:num_spikes, 
    plot(t_spikes(k),i_spikes(k),...    % plot spike rastergram
    '.k','Markersize',20); 
    hold on;
end;
set(gca,'Fontsize',16);
N_half=round(N/2);                      % beautify plot a bit
if N_half>1&N_half<N,
    set(gca,'Ytick',[1,N_half,N]);
else
    set(gca,'Ytick',[1:N]);
end;
    
axis([t_final-10,t_final,0,N+1]);        % show only last 5 time units
if max(max(delta))==min(min(delta)),
    delta_str=num2str(delta(1,2));
    title(['$\delta=$',delta_str,':'],'Fontsize',20);
end;
xlabel('$t$ [units of $T$]','Fontsize',20);
hold off;
shg;
                                        
toc     % Matlab computes execution time for statements bracketed by
        % "tic" and "toc". This can be useful when you play with the
        % program.
