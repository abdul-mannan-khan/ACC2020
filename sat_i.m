function [fx,dfx,ifx] = sat_i(u,m,M)
    if abs(u)<=m
        fx = u;
        dfx = 1;
        ifx = u^2/2;
    elseif u>m
        w = pi*(u-m) / (2*(M-m));
        fx = m+2*(M-m)*atan(w)/pi;
        dfx = (1+w^2)^(-1);
        ifx = m*(u-m/2)+(2*(M-m)/pi)^2*iatan(w);
    else
        w = pi*(u+m) / (2*(M-m));
        fx = -m+2*(M-m)*atan(w)/pi;
        dfx = (1+w^2)^(-1);
        ifx = -m*(u+m/2)+(2*(M-m)/pi)^2*iatan(w);
    end
end

function out = iatan(u)
    out = u*atan(u) - log(u^2 + 1)/2;
end