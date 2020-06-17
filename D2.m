function out = D2(xi)
    global delta
    p = xi(1:3);
    v = xi(4:6);
    R = reshape(xi(7:15),[3 3]);
    q = xi(16:18);
    t = xi(end);
    [pd,vd,ddot_pd] = reference(t);
    e = [p-pd;v-vd];
    x = R'*rd(ddot_pd,e);
    if V2(x,q)-minV2(x) >= delta
        out = 1;
    else
        out = 0;
    end
end