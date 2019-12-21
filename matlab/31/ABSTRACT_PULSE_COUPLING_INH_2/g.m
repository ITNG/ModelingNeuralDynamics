function g=g(phi);

a=4;
epsilon=0.4;

g=-epsilon*phi.*tanh((1-phi)*a)/tanh(a);