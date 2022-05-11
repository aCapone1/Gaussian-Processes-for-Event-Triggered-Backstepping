function dxdt = diff_olrt_gp_datagen(t,x,Xd1, Xd2, Xd3, Yd1, Yd2,Yd3,...
    postF1,postF2,postF3)

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    z11 = x(4);
    z12 = x(5);
    z21 = x(6);
    z22 = x(7);
    xi1 = x(8);
    xi2 = x(9);
    %Last compensation term is set to zero
    xi3 = 0;

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
    
    %Passive component
    F0 = 0;
    F1 = (-N*sin(x1)-B*x2)/D;
    F2 = (-Km*x2 - H*x3)/M;

    %Active/coupling component
    G0 = 1;
    G1 = 1/D;
    G2 = 1/M;

    %Desired signals and derivatives
    x1d = sin(2*pi*t);
    x1ddot = (2*pi)*cos(2*pi*t);
    x1dDot = -(2*pi)^2*sin(2*pi*t);
    x1dDdot = -(2*pi)^3*cos(2*pi*t);
    x2d = z11;
    x2ddot = omega_f*z12;
    x3d = z21;
    x3ddot = omega_f*z22;

    %System estimates
    [muF0 s1] = gp(hyp1, @infGaussLik, ...
     meanfunc, covfunc, likfunc, Xd1, postF1, x1);

    [muF1 s2] = gp(hyp2, @infGaussLik, ...
     meanfunc, covfunc, likfunc, Xd2, postF2, [x1,x2]);

    [muF2 s3] = gp(hyp3, @infGaussLik, ...
     meanfunc, covfunc, likfunc, Xd3, postF3, [x1,x2,x3]);

    %Pseudocontrol signals and derivatives
    alpha1 = 1/G0*(-K1*(x1-x1d) + x1ddot - muF0);
    alpha2 = 1/G1*(-K2*(x2-x2d) + x2ddot -muF1 -G0*((x1-x1d)-xi1));
    u = 1/G2*(-K3*(x3-x3d) + x3ddot - muF2 -G1*((x2-x2d)-xi2));

    
    xi1dot = -K1*xi1 + G0*(x2d - alpha1) +G0*xi2;
    xi2dot = -K2*xi2 + G1*(x3d - alpha2) +G1*xi3;
    z11dot = omega_f*z12;
    z12dot = -2*zeta*omega_f*z12 - omega_f*(z11 - alpha1);
    z21dot = omega_f*z22;
    z22dot = -2*zeta*omega_f*z22 - omega_f*(z21 - alpha2);

    dx1dt = F0 + G0*x2;
    dx2dt = F1 + G1*x3;
    dx3dt = F2 + G2*u;

    dxdt = [dx1dt;dx2dt;dx3dt;z11dot;z12dot;z21dot;z22dot;xi1dot;xi2dot]; 

end