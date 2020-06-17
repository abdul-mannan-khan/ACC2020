function out = rho(e)
    global rho_max B
    if e ~= 0
        [fx,dfx] = Vp(e);
        out = norm(B'*dfx')/sqrt(fx);%max([rho_max/3, norm(B'*dfx')/sqrt(fx)]);
    else
        out = rho_max;
    end