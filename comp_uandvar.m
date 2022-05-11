function [u, cov1, cov2, cov3] = comp_uandvar(t_sim,X_sim,Xd1,Yd1,Xd2,Yd2,Xd3,Yd3)

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
    omega_f=par(10);
    zeta=par(11);

    meanfunc = hyperpar.meanfunc;
    covfunc = hyperpar.covfunc;
    likfunc = hyperpar.likfunc;
    hyp1 = hyperpar.hyp1;
    hyp2 = hyperpar.hyp2;
    hyp3 = hyperpar.hyp3;
    
    u = zeros(size(t_sim));
    cov1 = zeros(size(t_sim));
    cov2 = zeros(size(t_sim));
    cov3 = zeros(size(t_sim));
    
    for i=1:length(t_sim)

        t = t_sim(i);
        x1 = X_sim(i,1);
        x2 = X_sim(i,2);
        x3 = X_sim(i,3);
        z11 = X_sim(i,4);
        z12 = X_sim(i,5);
        z21 = X_sim(i,6);
        z22 = X_sim(i,7);
        xi1 = X_sim(i,8);
        xi2 = X_sim(i,9);

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
        x2d = z11;
        x2ddot = omega_f*z12;
        x3d = z21;
        x3ddot = omega_f*z22;

%         %System estimates
        [muF0 s1] = gp(hyp1, @infGaussLik, ...
         meanfunc, covfunc, likfunc, Xd1, Yd1, x1);

        [muF1 s2] = gp(hyp2, @infGaussLik, ...
         meanfunc, covfunc, likfunc, Xd2, Yd2, [x1,x2]);

        [muF2 s3] = gp(hyp3, @infGaussLik, ...
         meanfunc, covfunc, likfunc, Xd3, Yd3, [x1,x2,x3]);

        %Pseudocontrol signals and derivatives
        u(i) = 1/G2*(-K3*(x3-x3d) + x3ddot - muF2 -G1*((x2-x2d)-xi2));
        cov1(i) = s1;
        cov2(i) = s2;
        cov3(i) = s3;
    end
end