function [Xdata_1_perm, Xdata_2_perm, Xdata_3_perm, ...
    Ydata_1_perm, Ydata_2_perm, Ydata_3_perm,...
    hyp1, hyp2, hyp3,postF1,postF2,postF3] = comp_hyperpar(npoints)


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
    std_noise=par(12);
    
    meanfunc = hyperpar.meanfunc;
    covfunc = hyperpar.covfunc;
    likfunc = hyperpar.likfunc;
    
    % hyp.cov = [ log(ell_1)
    %               log(ell_2)
    %               ...
    %               log(ell_D)
    %               log(sf) ]
    
    hyp1 = hyperpar.hyp1;
    hyp2 = hyperpar.hyp2;
    hyp3 = hyperpar.hyp3;


    Xdata_1_perm = zeros(npoints,1);
    Xdata_2_perm = zeros(npoints,2);
    Xdata_3_perm = zeros(npoints,3);

    %GENERATE DATA
    %GP PARAMETERS FOR SIMULATION
    Xd1 = rand(3,1);
    Xd2 = rand(3,2);
    Xd3 = rand(3,3);
    Yd1 = zeros(3,1);
    Yd2 = zeros(3,1);
    Yd3 = zeros(3,1);

    %Initial conditions
    [x10,x20,x30,z110,z120,z210,z220,xi10,xi20]...
        = get_initial_cond_datagen(Xd1, Xd2, Xd3, Yd1, Yd2,Yd3);

    Init_cond = [x10,x20,x30,z110,z120,z210,z220,xi10,xi20];

    [~, ~, postF1] = gp(hyp1, @infGaussLik, meanfunc, covfunc, ...
            likfunc, Xd1, Yd1);

    [~, ~, postF2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, ...
            likfunc, Xd2, Yd2);

    [~, ~, postF3] = gp(hyp3, @infGaussLik, meanfunc, covfunc, ...
            likfunc, Xd3, Yd3);

%     [t_data,x_data]=ode45(@(t,x) ...
%         OLRT_GPs_diff_CF(t,x,Xd1, Xd2, Xd3, Xd1, Xd2, Xd3, ...
%     Yd1, Yd2, Yd3, Yd1, Yd2, Yd3, postF1, postF2, postF3,...
%     postF1, postF2, postF3,  hyp1,hyp2,hyp3,omega_f,zeta), [0,0.5],Init_cond);

    [t_data,x_data]=ode45(@(t,x) ...
        diff_olrt_gp_datagen(t,x,Xd1, Xd2, Xd3, Yd1, Yd2,Yd3,...
    postF1,postF2,postF3), [0,1],Init_cond);

    %Initialize sigmamax
    sum_sigma_max = 0;

    Ydata_1 = zeros(npoints,1);
    Ydata_2 = zeros(npoints,1);
    Ydata_3 = zeros(npoints,1);

    perm = randperm(size(x_data,1));
    
    x_data_perm = x_data(perm,:);
    
    for i=1:npoints
        Xdata_1_perm(i,:) = x_data_perm(i,1);
        Ydata_1_perm(i,:) = 0 +randn*std_noise;
        Xdata_2_perm(i,:) = x_data_perm(i,1:2);
        Ydata_2_perm(i,:) = (-N*sin(x_data_perm(i,1))-B*x_data_perm(i,2))/D...
        +randn*std_noise;
        Xdata_3_perm(i,:) = x_data_perm(i,1:3);
        Ydata_3_perm(i,:) = (-Km*x_data_perm(i,2) - H*x_data_perm(i,3))/M...
        +randn*std_noise;
    end
    
%     % train hyperparameters. Normally commented out
%     [~, hyp1] = gp(hyp1, @infGaussLik, meanfunc, covfunc, ...
%         likfunc, Xdata_1_perm, Ydata_1_perm);
%     [~, hyp2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, ...
%         likfunc, Xdata_2_perm, Ydata_2_perm);
%     [~, hyp3] = gp(hyp3, @infGaussLik, meanfunc, covfunc, ...
%         likfunc, Xdata_3_perm, Ydata_3_perm);
    
    %Compute Post structures
    [~, ~, postF1] = gp(hyp1, @infGaussLik, meanfunc, covfunc, ...
        likfunc, Xdata_1_perm, Ydata_1_perm);
    [~, ~, postF2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, ...
        likfunc, Xdata_2_perm, Ydata_2_perm);
    [~, ~, postF3] = gp(hyp3, @infGaussLik, meanfunc, covfunc, ...
        likfunc, Xdata_3_perm, Ydata_3_perm);
  
end