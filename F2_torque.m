function [dxi,omegad] = F2_torque(xi,varargin)
    global r g m k2 motor_pole
    p = xi(1:3);
    v = xi(4:6);
    R = reshape(xi(7:15),[3 3]);
    q = xi(16:18);
    omega = xi(19:21);
    t = xi(22);
    %V = xi(20);
    
    [pd,dot_pd,ddot_pd,d3dot_pd] = reference(t);
    
    e = [p-pd;v-dot_pd];

    if nargin > 1 && varargin{1} == 1
        tau = varargin{2};
        [T,omegad] = openloop(tau);
    else
        T = T_2(ddot_pd,e);
        omegad = omega_2(d3dot_pd,ddot_pd,e,R,q);
    end
    dp = v;
    dv = R*r*T/m+g;
    dR = R*S(omega);
    dq = zeros(3,1);
    domega = 2*pi*motor_pole*(omegad-omega);
    dxi = [dp;dv;dR(:);dq;domega;1];
end



