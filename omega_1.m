function out = omega_1(d3dot_pd,ddot_pd,e,R,q)
    global r g alpha k1 m
    [f1,Deta] = eta(ddot_pd,e);
    out = -S(r)^2*(q*k1*S(r)*R'*rd(ddot_pd,e)+...
        q*nu(ddot_pd,e,R)/alpha+...
        S(R'*rd(ddot_pd,e))*R'*Deta*[d3dot_pd;e(4:6);R*r*T_1(ddot_pd,e,R)/m+g-ddot_pd]);
end

function out = nu(ddot_pd,e,R)
    global B r
    [~,gradVp] = Vp(e);
    out = -norm(eta(ddot_pd,e))*S(r)*R'*B'*gradVp';
end