function out = T_1(ddot_pd,e,R)
    global r m
    out = m*r'*R'*eta(ddot_pd,e);
end