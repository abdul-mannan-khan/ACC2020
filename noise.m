function [t,j,xi] = noise(s)
    global r
    params(s);
    TSPAN = [0 10];
    JSPAN = [0 10];
    rule  = 1;
    p0 = zeros(3,1);
    v0 = zeros(3,1);
    R0 = eye(3);
    q01= 1;
    q02 = -r;
    z10 = [p0;v0;R0(:);q01;zeros(3,1);0];
    z20 = [p0;v0;R0(:);q02;zeros(3,1);0];
    options = odeset('maxstep',0.1);
    [t,j,xi] = HyEQsolver(@F,@G,@C,@D,[z10;z20;zeros(4,1)], TSPAN,JSPAN,rule,options);
end

function dxi = F(xi)
    z1 = xi(1:20);
    z2 = xi(21:42);
    s1 = xi(43);
    tau1 = xi(44);
    s2 = xi(45);
    tau2 = xi(46);
    dxi= [F1_torque(z1,s1,tau1);F2_torque(z2,s2,tau2);0;1;0;1];
end

function out = C(xi)
    z1 = xi(1:20);
    z2 = xi(21:42);
    s1 = xi(43);
    tau1 = xi(44);
    s2 = xi(45);
    tau2 = xi(46);
    out = (s1 == 1 || C1(z1)) & (s2 == 1 || C2(z2));
end

function g = G(xi)
    z1 = xi(1:20);
    z2 = xi(21:42);
    s1 = xi(43);
    tau1 = xi(44);
    s2 = xi(45);
    tau2 = xi(46);
    if Ds(s1,tau1,z1(end))
        s1p = 1-s1;
        tau1p = 0;
    else
        s1p = s1;
        tau1p = tau1;
    end
    if Ds(s2,tau2,z2(end))
        s2p = 1-s2;
        tau2p = 0;
    else
        s2p = s2;
        tau2p = tau2;
    end
    if s1 == 0 && D1(z1)
        g1 = G1(z1);
    else
        g1 = z1;
    end
    if s2 == 0 && D2(z2)
        g2 = G2(z2);
    else
        g2 = z2;
    end
    
    g = [g1;g2;s1p;tau1p;s2p;tau2p];
end

function out = D(xi)
    z1 = xi(1:20);
    z2 = xi(21:42);
    s1 = xi(43);
    tau1 = xi(44);
    s2 = xi(45);
    tau2 = xi(46);
    out = (s1 == 0 && D1(z1)) | (s2 == 0 && D2(z2)) | Ds(s1,tau1,z1(end)) | Ds(s2,tau2,z2(end));
end

function out = Ds(s,tau,t)
    global max_w t0 n_flips
    if t >= t0 && t < t0+n_flips*2*pi/max_w/2
        u = 1;
    else
        u = 0;
    end
    if u == 1 && s == 0
        out = 1;
    elseif tau >= n_flips*2*pi/max_w && s == 1
        out = 1;
    else
        out = 0;
    end
end