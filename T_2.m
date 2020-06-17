function out = T_2(ddot_pd,e)
    global m
    out = m*norm(eta(ddot_pd,e));
end