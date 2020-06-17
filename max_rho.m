function out = max_rho()    
    global K 
    out = 0;
    for I = 1:3
        ki = [K(I,I);K(I,I+3)];
        Mi1= zeros(2);
        Mi1(1,1) = -ki(2);
        Mi2= [ki(2)^3/ki(1) ki(2)^2/2;ki(2)^2/2 ki(1)*ki(2)/2];
        tilde_Pi = [1-ki(2)^2/ki(1),-ki(2)/2;-ki(2)/2,-ki(1)];
        out = out+max(svd((Mi1+Mi2)*tilde_Pi^(1/2))); 
    end