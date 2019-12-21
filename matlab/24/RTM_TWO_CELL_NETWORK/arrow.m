function arrow=arrow(A,B,C,D,x,y,v,epsilon,width,col);

% input:  A,B,C,D = limits of plotting window. The plotting window
%         is [A,B] by [C,D]. Use axis('square') when using this program.
%         (x,y) = location of tip of arrow. 
%         v = vector in the direction of arrow;
%         epsilon = length of the two line segments that make up the arrow,
%         normalized to be between 0 and 1.
%         width=linewidth
%         col=a string like '-k' or '-r' to specify the color of
%              the arrow
% Typically, set "hold" to "on" before using this code. 
% The effect of this code will be to set "hold" to "on".

u=zeros(2,1);
u(1)=v(1)/(B-A);
u(2)=v(2)/(D-C);
u=u/norm(u)*epsilon;

u_right=[cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*u;
u_left =[cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)]*u;

v_right=zeros(2,1);
v_right(1)=u_right(1)*(B-A);
v_right(2)=u_right(2)*(D-C);

v_left=zeros(2,1);
v_left(1)=u_left(1)*(B-A);
v_left(2)=u_left(2)*(D-C);


h=plot([x,x-v_right(1)], ...
       [y,y-v_right(2)],col,'Linewidth',width);
set(h,'clipping','off');
hold on;
h=plot([x,x-v_left(1)], ...
     [y,y-v_left(2)],col,'Linewidth',width);
set(h,'clipping','off');


 