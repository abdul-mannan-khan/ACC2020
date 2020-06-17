function out = G2(xi)
    out = xi;
    p = xi(1:3);
    v = xi(4:6);
    R = reshape(xi(7:15),[3 3]);
    q = xi(16:18);
    t = xi(19);
    [pd,vd,ddot_pd] = reference(t);
    e = [p-pd;v-vd];
    x = R'*rd(ddot_pd,e);
    if D2(xi)
        out(16:18) = qminV2(x);
    end
end