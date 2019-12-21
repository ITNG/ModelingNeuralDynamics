clear; clf;

tic
%a=100/21; tau_n=100/sqrt(3); % parameters in the FN-equation
a=5; tau_n=60;              % parameters in the FN-equation
i_ext_vec=-6+[0:500]/500*12;% values of I to explore
i_c=-4.291556094395530;
i_ext_vec=1.01*i_c+[0:500]/500*(0.98*i_c-1.01*i_c);

subplot(111);
real_part_old=-1000;        % This initializes the piece of the code
                            % that determines where the fixed point loses
                            % stability, see below.

n_black=0;                  % stable nodes
n_red=0;                    % stable spirals
n_green=0;                  % unstable spirals
n_blue=0;                   % unstable nodes
n_magenta=0;                % saddles

t_final=10000;              % we find the amplitude of the limit cycle
                            % by solving from time 0 to time t_final ...
dt=0.01; dt05=dt/2;         % ... with time step dt
m_steps=round(t_final/dt);
v=zeros(m_steps+1,1); 
n=zeros(m_steps+1,1);

n_cycle=0;
for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);

    f=@(v) v-v^3/3-a*v+i_ext;   % zeros of this function are fixed points
                                % of the FN system

    v_left=-1; v_right=1;
    while f(v_left)<0,
        v_left=v_left-1;
    end;                        % find v_left with f(v_left)>0
    while f(v_right)>0,
        v_right=v_right+1;
    end;                        % find v_right with f(v_right)<0

    while v_right-v_left>10^(-10),
        v_c=(v_left+v_right)/2;
        if f(v_c) >=0, 
            v_left=v_c;         % find a zero of f using bisection
        else
            v_right=v_c;
        end;
    end;

    v_c=(v_left+v_right)/2;     % fixed point of the FN system
    n_c=a*v_c;

    J=zeros(2,2);
    J(1,1)=1-v_c^2;
    J(1,2)=-a;
    J(2,1)=a/tau_n;
    J(2,2)=-1/tau_n;            % J = Jacobi matrix of the FN system
                                %     at the fixed point
    
    E=eig(J);
    real_part=real(E(1));
    if real_part>0 && real_part_old<0,
        i_c=(i_ext_vec(ijk)*(-real_part_old)+i_ext_vec(ijk-1)*real_part)/ ...
        (real_part-real_part_old)
    end;                        % This finds the value at which 
                                % the fixed point loses stability.
    real_part_old=real_part;
    
    
    if abs(imag(E(1)))>10^(-4), % non-real eigenvalues
        if real(E(1))<0,        % stable spiral
            n_red=n_red+1;      
            i_red(n_red)=i_ext;
            v_c_red(n_red)=v_c;
        end;
        if real(E(1))>0,        % unstable spiral
            n_green=n_green+1;
            i_green(n_green)=i_ext;
            v_c_green(n_green)=v_c;
        end;
    else
        if E(1)<0 & E(2)<0,     % stable node
            n_black=n_black+1;
            i_black(n_black)=i_ext;
            v_c_black(n_black)=v_c;
        end;
        if E(1)>0 & E(2)>0,     % unstable node
            n_blue=n_blue+1;
            i_blue(n_blue)=i_ext;
            v_c_blue(n_blue)=v_c;
        end;
        if E(1)>0 & E(2)<0,     % saddle
            n_magenta=n_magenta+1;
            i_magenta(n_magenta)=i_ext;
            v_c_magenta(n_magenta)=v_c;
        end;
    end;
    
    if abs(imag(E(1)))>10^(-4) && real(E(1))>0,     % unstable spiral
        v(1)=v_c; n(1)=n_c*1.05;
        for k=1:m_steps,
            v_inc=v(k)-v(k)^3/3-n(k)+i_ext;
            n_inc=(a*v(k)-n(k))/tau_n;
            v_tmp=v(k)+dt05*v_inc;
            n_tmp=n(k)+dt05*n_inc;
            v_inc=v_tmp-v_tmp^3/3-n_tmp+i_ext;
            n_inc=(a*v_tmp-n_tmp)/tau_n;
            v(k+1)=v(k)+dt*v_inc;
            n(k+1)=n(k)+dt*n_inc;
        end;
        vv=v(4*m_steps/5:m_steps+1); 
        n_cycle=n_cycle+1;
        i_cycle(n_cycle)=i_ext;
        max_cycle(n_cycle)=max(vv);
        min_cycle(n_cycle)=min(vv);
    end;
end;
if n_blue>0,
    plot(i_blue,v_c_blue,'-.b','Linewidth',4);   % unstable nodes
    hold on;
end;
if n_green>0,
    plot(i_green,v_c_green,'--g','Linewidth',4); % unstable spirals
    hold on;
end;
if n_red>0,
    plot(i_red,v_c_red,'.r','Markersize',20);     % stable spirals
    hold on;
end;
if n_black>0,
    plot(i_black,v_c_black,'-k','Linewidth',4); % stable nodes
    hold on;
end;
if n_magenta>0,
    plot(i_magenta,v_c_magenta,'--m','Linewidth',4);    % saddles
    hold on;
end;
if n_cycle>0,
    plot(i_cycle,max_cycle,':k','Linewidth',2);
    hold on;
    plot(i_cycle,min_cycle,':k','Linewidth',2);
    hold off;
end;

set(gca,'Fontsize',24);
xlabel('$I$', 'Fontsize',32);
ylabel('$v_\ast$', 'Fontsize',32);
set(gca,'Xtick',[-4.3, -4.26,-4.22]);
A=i_ext_vec(1); 
B=i_ext_vec(length(i_ext_vec));
axis([A,B,-2.5,2]); axis('square');
hold off;
shg;
toc
