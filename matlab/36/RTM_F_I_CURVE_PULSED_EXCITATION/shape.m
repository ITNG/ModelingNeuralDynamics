function shape=shape(t)

global alpha Period ave;

shape=(exp(alpha*cos(pi*t/Period).^2)-1)/ave;