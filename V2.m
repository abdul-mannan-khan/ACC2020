function [fx,dfx] = V2(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [Vxy,gradVy] = V(x,y)
%  x,y - unitary vectors
%  Vxy = V(x,y) = (1-r' x)/(1-r' x+k(1-y' x))
%  gradVy = dV/dx
% globals
%  kr - positive scalar
%  r - reference unitary vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global kr r
    fx = (1-r'*x)/(1-r'*x+kr*(1-y'*x));
    dfx = (kr*fx*y-(1-fx)*r)/(1-r'*x+kr*(1-y'*x));