function g=g(phi);
epsilon=0.75;
g=2*(1-phi)-(1+epsilon-sqrt((1+epsilon)^2-4*epsilon*(1-phi)))/epsilon;
