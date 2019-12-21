clear;  clf;

tau_m_0=10; J_0=0.02; g_0=0.15; tau_I_0=9; % baseline point

J_vec=J_0*0.8+[0:100]/100*(J_0*1.2-J_0*0.8);
tau_I_vec=tau_I_0*0.8+[0:100]/100*tau_I_0*0.4;
g_vec=g_0*0.8+[0:100]/100*g_0*0.4;


for k=1:length(tau_I_vec),
    P_vec(k)=P(tau_m_0,J_0,g_0,tau_I_vec(k)); S_vec(k)=S(tau_m_0,J_0,g_0,tau_I_vec(k));
end;

subplot(331);
plot(tau_I_vec,P_vec,'-k','Linewidth',2);
set(gca,'Fontsize',12);
ylabel('$P$','Fontsize',16);
title('varying $\tau_I$','Fontsize',16);
hold on;
plot(tau_I_0,P(tau_m_0,J_0,g_0,tau_I_0),'.r','Markersize',20);
hold off;
axis([tau_I_0*0.8,tau_I_0*1.2,20,40]);


subplot(334);
plot(tau_I_vec,S_vec,'-k','Linewidth',2);
set(gca,'Fontsize',12);
xlabel('$\tau_I$','Fontsize',16);
ylabel('$S$','Fontsize',16);
axis([tau_I_0*0.8,tau_I_0*1.2,0,0.1]);
hold on;
plot(tau_I_0,S(tau_m_0,J_0,g_0,tau_I_0),'.r','Markersize',20);
hold off;

shg;



for k=1:length(g_vec),
    P_vec(k)=P(tau_m_0,J_0,g_vec(k),tau_I_0);
    S_vec(k)=S(tau_m_0,J_0,g_vec(k),tau_I_0);
end;

subplot(332);
plot(g_vec,P_vec,'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[0.13:0.02:0.17]);
axis([0.8*g_0,1.2*g_0,20,40]);
title('varying $\overline{g}_{\rm syn}$','Fontsize',16);
hold on;
plot(g_0,P(tau_m_0,J_0,g_0,tau_I_0),'.r','Markersize',20);
hold off;


subplot(335);
plot(g_vec,S_vec,'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[0.13:0.02:0.17]);
xlabel('$\overline{g}_{\rm syn}$','Fontsize',16);
axis([0.8*g_0,1.2*g_0,0,0.1]);
hold on;
plot(g_0,S(tau_m_0,J_0,g_0,tau_I_0),'.r','Markersize',20);
hold off;
shg;


for k=1:length(J_vec),
    P_vec(k)=P(tau_m_0,J_vec(k),g_0,tau_I_0);
    S_vec(k)=S(tau_m_0,J_vec(k),g_0,tau_I_0);
end;

subplot(333);
plot(J_vec,P_vec,'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[0.017,0.02,0.023]);
axis([0.8*J_0,1.2*J_0,20,40]);
title('varying $J$','Fontsize',16);
hold on;
plot(J_0,P(tau_m_0,J_0,g_0,tau_I_0),'.r','Markersize',20);
hold off;


subplot(336);
plot(J_vec,S_vec,'-k','Linewidth',2);
set(gca,'Fontsize',12);
set(gca,'Xtick',[0.017,0.02,0.023]);
xlabel('$J$','Fontsize',16);
axis([0.8*J_0,1.2*J_0,0,0.1]);
hold on;
plot(J_0,S(tau_m_0,J_0,g_0,tau_I_0),'.r','Markersize',20);
hold off;
shg;

