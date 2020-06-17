function [T,omega] = openloop(tau)
    global max_w m g
    T = m*norm(g);
    if tau < 2*pi/max_w
        omega = [max_w;0;0];
    else
        omega = zeros(3,1);
    end