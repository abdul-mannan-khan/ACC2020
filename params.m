function params(varargin)
    global K m r g ell M Pi Mi1 Mi2 ee k1 B alpha rho_max k2 gg kr c1 c2 c3 delta max_w motor_pole A n_flips t0
    t0 = 1;
    max_w = 18.5; 
    motor_pole = 6;
    
    Qlqr = diag([ones(1,3) 0.01*ones(1,3)]);
    Rlqr = eye(3);
    A = [zeros(3),eye(3);zeros(3), zeros(3)];
    B = [zeros(3);eye(3)];
    K = -lqr(A,B,Qlqr,Rlqr);
    Pi= zeros(2,2,3);
    Mi1 = zeros(2,2,3);
    Mi2 = zeros(2,2,3);
    for I = 1:3
        k1 = K(I,I);
        k2 = K(I,I+3);
        Mi1(:,:,I) = diag([-k2/2 0]);
        Mi2(:,:,I) = [k1^3/(2*k1), k2^2/2;
                      k2^2/2,      k1*k2/2];
        Pi(:,:,I) = [-k2^2/(2*k1), -k2/2;
                    -k2/2, -k1];
    end
    if nargin>0
        alpha = varargin{1};
    else
        alpha = 1e3;
    end
    
    k1 = 1;
    r = [0;0;-1];
    g = [0;0;9.81];
    ell = ones(3,1)*0.19*norm(g);
    m = 0.216;
    M = 2*ell;
    ee = ones(3,1);
    rho_max = max_rho();
    k2 = 1;
    gg = 0;
    kr = 1;
    c1 = 1/(2*(1+kr+sqrt(1+kr*gg+kr^2)));
    c2 = 1/(2*(1+kr-sqrt(1+kr*gg+kr^2)));
    delta = 0.5*(1+gg)/(2/kr+(1+gg));
    Vstar = 2/(2+kr*(1+gg))+delta;
    c3 = (2*kr*(1-Vstar)*(1-gg))/(1+kr+sqrt(1+2*kr*gg+kr^2));
    n_flips = 0.5;
    