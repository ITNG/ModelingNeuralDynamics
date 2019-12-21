clear; clf;

c=1;
g_na=20;
g_k=10; 
g_k_slow=5;
g_l=8;
v_na=60;
v_k=-90;
v_l=-80;
tau_n=0.15;
tau_n_slow=20;

i_ext=7;

t_final=1500; dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); v=z; m=z; n=z; n_slow=z;
v(1)=-70; m(1)=m_inf(v(1)); n(1)=0.6; n_slow(1)=0;

num_cyc=0; % num_cyc counts how many times the trajectory has gone around
           % the limit cycle; see below for how that is done. 
for k=1:m_steps,
    
    v_inc=(g_na*m(k)*(v_na-v(k))+ ...
        g_k*n(k)*(v_k-v(k))+ ...
        g_k_slow*n_slow(k)*(v_k-v(k))+ ...
        g_l*(v_l-v(k))+i_ext)/c;
    n_inc=(n_inf(v(k))-n(k))/tau_n;
    n_slow_inc=(n_slow_inf(v(k))-n_slow(k))/tau_n_slow;
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    n_tmp=n(k)+dt05*n_inc;
    n_slow_tmp=n_slow(k)+dt05*n_slow_inc;
    
    v_inc=(g_na*m_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp*(v_k-v_tmp)+ ...
        g_k_slow*n_slow_tmp*(v_k-v_tmp)+ ...
        g_l*(v_l-v_tmp)+i_ext)/c;
    n_inc=(n_inf(v_tmp)-n_tmp)/tau_n;
    n_slow_inc=(n_slow_inf(v_tmp)-n_slow_tmp)/tau_n_slow;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    n(k+1)=n(k)+dt*n_inc;
    n_slow(k+1)=n_slow(k)+dt*n_slow_inc;
    if n_slow(k)>0.02 & n_slow(k+1)<=0.02,
        num_cyc=num_cyc+1;  % We say that the trajectory has gone
                            % around the cycle once when n_slow drops
                            % below 0.02. (That 0.02 is a good value
                            % to use here can of course not be known
                            % a priori, I figured that out by looking at
                            % the computed solution.) 
        k_vec(num_cyc)=k;   % Remember what k was when the present round 
                            % was completed. 
    end;
    
end;

subplot(111);
t=[0:m_steps]*dt;
k1=k_vec(num_cyc-1); k2=k_vec(num_cyc); % We plot only the last passage
                                        % through the limit cycle, so that
                                        % we'll see the cycle, but not the
                                        % approach to it, in the figure.
ind=k1:k2;
v=v(ind); n=n(ind); n_slow=n_slow(ind);
plot3(v,n,n_slow,'-k','Linewidth',2);
set(gca,'Fontsize',24);
xlabel('$v$','Fontsize',32);
ylabel('$n$','Fontsize',32);
zlabel('$n_{\rm slow}$','Fontsize',32);
axis([-70,0,0,0.7,0,0.06]);
axis('square');

% This code does not add the little arrow that you see in the picture, 
% indicating the direction of motion. I found it too much of a pain
% to add that using Matlab code --- I just added it using "Insert Arrow" 
% in Matlab. 

shg;