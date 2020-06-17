function out = rd(ddot_pd,e)
    feta = eta(ddot_pd,e);
    out  = feta/norm(feta);
end