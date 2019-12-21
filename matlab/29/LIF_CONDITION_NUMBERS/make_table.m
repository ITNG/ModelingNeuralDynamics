tau_m=10;
J=0.02;
g=0.15;
tau_I=9;

(P(tau_m,0.99*J,g,tau_I)-P(tau_m,J,g,tau_I))/P(tau_m,J,g,tau_I)*100
(P(tau_m,J,g*1.01,tau_I)-P(tau_m,J,g,tau_I))/P(tau_m,J,g,tau_I)*100
(P(tau_m,J,g,tau_I*1.01)-P(tau_m,J,g,tau_I))/P(tau_m,J,g,tau_I)*100


g=2;
tau_I=1;

(P(tau_m,0.99*J,g,tau_I)-P(tau_m,J,g,tau_I))/P(tau_m,J,g,tau_I)*100
(P(tau_m,J,g*1.01,tau_I)-P(tau_m,J,g,tau_I))/P(tau_m,J,g,tau_I)*100
(P(tau_m,J,g,tau_I*1.01)-P(tau_m,J,g,tau_I))/P(tau_m,J,g,tau_I)*100