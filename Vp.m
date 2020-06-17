function [fx,dfx] = Vp(e)
   global K Pi ell M ee;
    fx = 0;
    dfx =  zeros(1, 6);
    for I = 1:3
        k = [K(I,I);K(I,I+3)];
        ei = [e(I);e(I+3)];
        [a,da] = Vpi(ei,Pi(:,:,I),k,ell(I),M(I),ee(I));
        fx = fx+a;
        dfx(I) = da(1);
        dfx(I+3) = da(2);
    end