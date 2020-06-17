function out = qminV2(x)
    global r gg
    if all(x == r)
        out = rand(3,1)*2-1;
        out = oout/norm(out);
    elseif r'*x >= -gg
        out = -x;
    elseif r'*x > -1 && r'*x<-gg
        out = sigma(r'*x)*PT(x)*r/norm(PT(x)*r)+alpha(r'*x)*x;
    else
        v = rand(3,1)*2-1;
        out = sqrt(1-gg^2)*PT(r)*v/norm(PT(r)*v)+gg*r;
    end
end

function out = PT(x)
    out = eye(3)-x*x';
end

function out = alpha(v)
    global gg
    out = gg*v-sqrt((1-v^2)*(1-gg^2));
end

function out = sigma(v)
    global gg
    out = gg*sqrt(1-v^2)+v*sqrt(1-gg^2);
end