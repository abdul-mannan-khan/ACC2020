function [fx,dfx] = sigma(v)
    global ell M
    N = numel(v);
    fx = zeros(N,1);
    dfx = zeros(N);
    for I = 1:numel(v)
        [f1,df1] = sat_i(v(I),ell(I),M(I));
        fx(I) = f1;
        dfx(I,I) = df1;
    end
    
end