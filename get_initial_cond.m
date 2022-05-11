function [z110,z120,z210,z220,xi10,xi20] = ...
    init(x10,x20,x30, Xd1, Xd2, Xd3, Yd1, Yd2,Yd3)

%     x10=0.1;
%     x20=6.28;
%     x30=0;

%     x10=0;
%     x20=0;
%     x30=0;


    %System parameters
    [par, hyperpar] = get_parameters();
    M = par(1);
    N = par(2);
    B = par(3);
    D = par(4);
    Km = par(5);
    H = par(6);
    K1 = par(7);
    K2 = par(8);
    K3 = par(9);
    omega_f = par(10);
    zeta = par(11);

    meanfunc = hyperpar.meanfunc;
    covfunc = hyperpar.covfunc;
    likfunc = hyperpar.likfunc;
    hyp1 = hyperpar.hyp1;
    hyp2 = hyperpar.hyp2;
    hyp3 = hyperpar.hyp3;
    
    %Active/coupling component
    G0 = 1;
    G1 = 1/D;
    G2 = 1/M;
    
    
    %System estimates
    muF00 = gp(hyp1, @infGaussLik, ...
    meanfunc, covfunc, likfunc, Xd1, Yd1, x10);

    muF10 = gp(hyp2, @infGaussLik, ...
    meanfunc, covfunc, likfunc, Xd2, Yd2, [x10,x20]);


    %Starting conditions of signals
    t0 = 0;
    
    x1d0 = sin(2*pi*t0);
    x1ddot0 = (2*pi)*cos(2*pi*t0);
    
    % xi here refers to the greek letter xi. Its initial condition is set
    % to zero, as detailed in the papers by Capone and Farrell
    xi10 = 0;
    xi20 = 0;
    z110 = 1/G0*(-K1*(x10-x1d0) + x1ddot0 - muF00);
    x2d0 = z110;
    z120 = 0;
    x2ddot0 = omega_f*z120;
    z210 = 1/G1*(-K2*(x20-x2d0) + x2ddot0 -muF10 -G0*((x10-x1d0)-xi10));
    z220 = 0;

end