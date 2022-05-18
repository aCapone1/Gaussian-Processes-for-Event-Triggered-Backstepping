clear
close all
clc


% Add folders to path
addpath('gpml','visualization','results')

% Initialize GPML toolbox
startup

% Number of repetitions
nreps = 100;

%Number of training points per dimension (Number of training 
%points = n_points^3)
n_points_untrig = 400;
n_points_trig = 1;

for rep = 1:nreps
    
    % Seeding random variables for reproducibility
    rng(rep)
    
    %Starting and ending time of simulation
    TStart = 0;
    TEnd = 10; % 10
    
    %% Simulation without online learning
    %Compute hyperparameters for first set of training points
    [Xd1, Xd2, Xd3, Yd1, Yd2, Yd3, hyp1, hyp2, hyp3,postF1,postF2,postF3]...
    = comp_hyperpar(n_points_untrig);


    % Initial conditions
    x10= 10*randn; %4.9570; 
    x20= 10*randn; %28.2840; 
    x30= 10*randn; %5.5908; 
    
    % Compute initial conditions for command filter
    [z110,z120,z210,z220,xi10,xi20]...
        = get_initial_cond(x10,x20,x30, Xd1, Xd2, Xd3, Yd1, Yd2,Yd3);
    Init_cond = [x10,x20,x30,z110,z120,z210,z220,xi10,xi20];

    % simulate without online learning
    [time_untrig,state_untrig]=ode45(@(t,x) ...
                diff_olrt_gp(t,x,Xd1, Xd2, Xd3,... % Yd1, Yd2,Yd3,...
            postF1,postF2,postF3), [TStart,TEnd],Init_cond);

%     %Compute u for untriggered setting
%     [u_untrig, ~, ~, ~] = comp_uandvar(t_untrig,x_untrig,Xd1,Yd1,Xd2,Yd2,Xd3,Yd3,...
%         hyp1,hyp1,hyp3);

    u_untrig = [];

    %% Simulation with event-triggered online learning

    %Compute hyperparameters for first set of training points
    [Xd1, Xd2, Xd3, Yd1, Yd2, Yd3, hyp1, hyp2, hyp3,postF1,postF2,postF3]...
    = comp_hyperpar(n_points_trig);

    x_data = Xd3;


    %Initial conditions
    [z110,z120,z210,z220,xi10,xi20]...
        = get_initial_cond(x10 ,x20 ,x30 , Xd1, Xd2, Xd3, Yd1, Yd2,Yd3);
    Init_cond = [x10,x20,x30,z110,z120,z210,z220,xi10,xi20];

    time = [];
    triggertime = [];
    state = [];
    input = [];
    est_err = [];
    
    % Flag is set to 1 if GP update was successful, 0 otherwise (if a data
    % point could not be added due to poor matrix condition)
    while isempty(time) || time(end) ~= TEnd

        Opt = odeset('Events', @(t,x) trigger(t,x,Xd1, Xd2, Xd3, Yd1, Yd2,Yd3,...
            postF1,postF2,postF3,gp_update_flag));
        
        if isempty(triggertime)
            t = TStart;
            x = Init_cond;
        else
            [t,x]=ode45(@(t,x) ...
                    diff_olrt_gp(t,x,Xd1, Xd2, Xd3, ...%Yd1, Yd2,Yd3,...
                postF1,postF2,postF3), [TStart,TEnd],Init_cond, Opt);
        end
        x_data_new = x(end,1:3);

        %Compute u
        [u, cov1, cov2, cov3] = comp_uandvar(t,x,Xd1,Yd1,Xd2,Yd2,Xd3,Yd3);

        [Xd1, Xd2, Xd3, Yd1, Yd2, Yd3, postF1,postF2,postF3,gp_update_flag]...
            = update_gp(x_data_new,Xd1, Xd2, Xd3, Yd1, Yd2, Yd3);

        TStart = t(end);
        Init_cond = x(end,:);
        % concatenate simulated segments of state and input
        if TStart == TEnd
            time = [time; t];
            state = [state; x];
            input = [input; u];
            est_err = [est_err; sqrt(cov1+cov2+cov3)];
        else
            time = [time; t(1:end-1)];
            state = [state; x(1:end-1,:)];
            input = [input; u(1:end-1,:)];
            est_err = [est_err; ...
                sqrt(cov1(1:end-1)+cov2(1:end-1)+cov3(1:end-1))];
            triggertime = [triggertime; t(end)];
        end

    end

    save("results/results_evnttrig" + int2str(rep))
%     disp('Plotting results...')
%     plot_all_results(time,state,time_untrig,state_untrig,input,est_err, triggertime)
%     disp('Press any key to continue')
%     pause;
    close all
end
    
%% Compute mean and variance of normalized L2 error
[l2err_trig, l2err_untrig, Ndata_trig] = compute_l2err(nreps);

% Normalize L2 errors
norm_factor = max(max(l2err_trig, l2err_untrig));
l2err_trig = 1/norm_factor*l2err_trig;
l2err_untrig = 1/norm_factor*l2err_untrig;

% compute mean L2 errors
avg_L2trig = mean(l2err_trig);
avg_L2untrig = mean(l2err_untrig);

plot_eventtrig_res(time,state,input, est_err,'triggered',triggertime)
plot_eventtrig_res(time_untrig,state_untrig,u_untrig, est_err,'untriggered',[])

%% Triggering function
function [trigger_val, isterminal, direction] = trigger(t, x, Xd1, Xd2, Xd3, ...
    Yd1, Yd2,Yd3,postF1,postF2,postF3,gp_update_flag)


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

    %System parameters and hyperparameters
    [par, hyperpar] = get_parameters();
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
    
    %determine smallest gain
    Kmin = min(K1,min(K2,K3));
    
    %Desired signals and derivatives
    x1d = sin(2*pi*t);
    x1ddot = (2*pi)*cos(2*pi*t);
    %x1dDot = -(2*pi)^2*sin(2*pi*t);
    %x1dDdot = -(2*pi)^3*cos(2*pi*t);
    x2d = z11;
    x2ddot = omega_f*z12;
    x3d = z21;
    x3ddot = omega_f*z22;
    
    beta = 2; %2
    err = sqrt((x1-x1d)^2 + (x2-x2d)^2 + (x3-x3d)^2);
    scaling_factor = 2*sqrt(beta)*1/Kmin;
   

    %System estimates
    [~, s1] = gp(hyp1, @infGaussLik, ...
     meanfunc, covfunc, likfunc, Xd1, postF1, x1);

    [~, s2] = gp(hyp2, @infGaussLik, ...
     meanfunc, covfunc, likfunc, Xd2, postF2, [x1,x2]);

    [~, s3] = gp(hyp3, @infGaussLik, ...
     meanfunc, covfunc, likfunc, Xd3, postF3, [x1,x2,x3]);

    sig_max = (exp(hyp1.lik) + exp(hyp2.lik) + exp(hyp3.lik));
    norm_var = sqrt(s1 + s2 + s3);
%     trigger_val = norm_var > 1.1*sqrt(sig_max);
    if gp_update_flag == 0
        add_factor = -1e-5;
    else
        add_factor = 0;
    end
    
    if err > scaling_factor*sqrt(sig_max)
        trigger_val = scaling_factor*norm_var - err + add_factor;
    else
        trigger_val = err - scaling_factor*sqrt(sig_max);
    end
    isterminal = 1;   % Stop the integration
    direction  = 1;
end
