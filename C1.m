function out = C1(xi)
    global delta r
    p = xi(1:3);
    v = xi(4:6);
    R = reshape(xi(7:15),[3 3]);
    q = xi(16);
    t = xi(17);
    [pd,vd,ddot_pd] = reference(t);
    e = [p-pd;v-vd];
    x = R'*rd(ddot_pd,e);
    if -2*q*r'*x <= delta
        out = 1;
    else
        out = 0;
    end
end