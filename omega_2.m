function out = omega_2(d3dot_pd,ddot_pd,e,R,q)
    global r g k2 m
    [f1,Deta] = eta(ddot_pd,e);
    dot_eta = Deta*[d3dot_pd;e(4:6);R*r*T_2(ddot_pd,e)/m+g-ddot_pd];
    [~,gradV2] = V2(R'*rd(ddot_pd,e),q);
    out = S(R'*rd(ddot_pd,e))*(R'*dot_eta/norm(f1)...
        +(k2+nu(ddot_pd,e))*gradV2);
end

function out = nu(ddot_pd,e)
    global alpha c1 c3
    out = norm(eta(ddot_pd,e))*rho(e)/(alpha*sqrt(c3*c1));
end