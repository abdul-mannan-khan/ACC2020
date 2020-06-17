function [fx, dfx] = Vpi(e,P,k,m,M,ee)
    u = k'*e;
    [s,ds,is] = sat_i(u,m,M);
    v = [s;e(2)];
    fx = 0.5*ee*v'*P*v+ee*is;
    dfx = ee*([s,e(2)]*P*[ds*k';0 1]+s*k');
end
