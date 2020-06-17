function [fx,dfx] = eta(ddot_pd,e)
    global K g
    [s,ds] = sigma(K*e);
    fx = s-g+ddot_pd;
    dfx = [eye(3),ds*K];
end