function [Xdata_1, Xdata_2, Xdata_3, Ydata_1, Ydata_2, Ydata_3,...
    postF1,postF2,postF3,gp_update_flag] = update_gp(x_data_new,...
    Xdata_1, Xdata_2, Xdata_3, Ydata_1, Ydata_2, Ydata_3)
    
    npoints = size(x_data_new,1);
    gp_update_flag = 1; % default
    
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
    std_noise=par(12);
    
    meanfunc = hyperpar.meanfunc;
    covfunc = hyperpar.covfunc;
    likfunc = hyperpar.likfunc;
    hyp1 = hyperpar.hyp1;
    hyp2 = hyperpar.hyp2;
    hyp3 = hyperpar.hyp3;
    

    for j=1:size(x_data_new,1)


            Xdata_1(end+j,:) = x_data_new(j,1);
            Ydata_1(end+j,:) = 0 +(0.5-rand)*std_noise;
            Xdata_2(end+j,:) = x_data_new(j,1:2);
            Ydata_2(end+j,:) = (-N*sin(x_data_new(j,1))-B*x_data_new(j,2))/D...
                +(0.5-rand)*std_noise;
            Xdata_3(end+j,:) = x_data_new(j,1:3);
            Ydata_3(end+j,:) = (-Km*x_data_new(j,2) - H*x_data_new(j,3))/M...
                +(0.5-rand)*std_noise;
    end
    %Compute Post structures
    try
        [~, ~, postF1] = gp(hyp1, @infGaussLik, meanfunc, covfunc, ...
            likfunc, Xdata_1, Ydata_1);
    catch
        Xdata_1 = Xdata_1(1:end-1,:);
        Ydata_1 = Ydata_1(1:end-1,:);
        
        [~, ~, postF1] = gp(hyp1, @infGaussLik, meanfunc, covfunc, ...
            likfunc, Xdata_1, Ydata_1);
        gp_update_flag = 0;
    end

    try
        [~, ~, postF2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, ...
        likfunc, Xdata_2, Ydata_2);
    catch
        Xdata_2 = Xdata_2(1:end-1,:);
        Ydata_2 = Ydata_2(1:end-1,:);
        
        [~, ~, postF2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, ...
        likfunc, Xdata_2, Ydata_2);
        gp_update_flag = 0;
    end

    try
        [~, ~, postF3] = gp(hyp3, @infGaussLik, meanfunc, covfunc, ...
        likfunc, Xdata_3, Ydata_3);
    catch
        Xdata_3 = Xdata_3(1:end-1,:);
        Ydata_3 = Ydata_3(1:end-1,:);
        
        [~, ~, postF3] = gp(hyp3, @infGaussLik, meanfunc, covfunc, ...
    likfunc, Xdata_3, Ydata_3);
        gp_update_flag = 0;
    end

end