function out = minV2(x)
    global r gg kr
    if r'*x >= -gg
        out = (1-r'*x)/(1-r'*x+2*kr);
    else
        out = (1-r'*x)/(1-r'*x+kr*(1-alpha(r'*x)));
    end
end

function out = alpha(v)
    global gg
    out = gg*v-sqrt((1-v^2)*(1-gg^2));
end