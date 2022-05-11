function [] = plot_all_results(t,x,t_untrig,x_untrig,u,est_err, triggertime)
%     close all

    par = get_parameters();
    K1 = par(7);
    K2 = par(8);
    K3 = par(9);
    omega_f=par(10);
    
    %Desired signals and derivatives
    x1d = sin(2*pi*t);
    z11 = x(:,4);
    z21 = x(:,6);
    x2d = z11;
    x3d = z21;
    
    err = ((x(:,1)-x1d).^2 + (x(:,2)-x2d).^2 + (x(:,3)-x3d).^2).^(0.5);
    %determine smallest gain
    Kmin = min(K1,min(K2,K3));
    beta = 2;
    scaling_factor = Kmin/(2*sqrt(beta));

    
    TEnd = t(end);
    err1 = abs(x(:,1)-sin(2*pi*t));
    err2 = abs(x(:,2)-2*pi*cos(2*pi*t));
    err3 = abs(x(:,3)-x(:,6));
    errtot = err1+err2+err3;
    
    err1_untrig = abs(x_untrig(:,1)-sin(2*pi*t_untrig));
    err2_untrig = abs(x_untrig(:,2)-2*pi*cos(2*pi*t_untrig));
    err3_untrig = abs(x_untrig(:,3)-x_untrig(:,6));
    errtot_untrig = err1_untrig+err2_untrig+err3_untrig;
    
    
    figure('Units','inches','Position',[0 0 10 4.5],...
    'PaperPositionMode','auto');
    set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
    'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
    'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
    'latex')
    hold on
    plot(t,errtot, 'b-', 'LineWidth', 2);
    hold on
    plot(t_untrig,errtot_untrig, 'r-', 'LineWidth', 2);
    %hold on
    %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
    grid off
    set(gca,'YScale','log')
    axis([0 TEnd 0.05 1000]) %-0.1*max(errtot) ...
%     1.1*max(errtot)])
%     legend({'Adaptive bs.'
%     'Bs. with GPs, $$N=50$$'},...% '$$N=4500$$'},...
%     'FontName','CMU Serif');
%     title({'$$x-x_d$$'}, 'FontSize', 22)
    xlabel('Time t (sec)', 'FontSize', 22)
    ylabel('Tracking err. $$\Vert e\Vert$$', 'FontSize', 22)
    set(gca,'fontsize',22)
    set(gca,'XTick', linspace(0,TEnd,6))
    scale = 0.1;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
    print('results/TotalError','-depsc2')
    
    
%     %NO INITIALIZATION
%     N3N5=CalcL2Norm(t1(470:end),x1(470:end,3)-x1(470:end,6));
%     N3N50=CalcL2Norm(t1(511:end),x1(511:end,3)-x1(511:end,6));
%     N3Ad=CalcL2Norm(Error3Kwan(59000:end,1),Error3Kwan(59000:end,2));
    
    if ~isempty(u)
        figure('Units','inches','Position',[0 0 9 4.5],...
        'PaperPositionMode','auto');
        set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
            'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
            'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
            'latex')
        hold on
        plot(t,u, 'b', 'LineWidth', 2);
        %hold on
        %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
        grid off
        axis([0 TEnd -1.4*abs(min(u)) ...
            1.6*abs(max(u))])
    %     legend({'Adaptive bs.','$$N=5$$',...
    %         '$$N=50$$'},...% '$$N=4500$$'},...
    %         'FontName','CMU Serif');
    %     title({'Simulation Results'}, 'FontSize', 22)
        xlabel('Time t (sec)', 'FontSize', 22)
        ylabel('Input $$u$$', 'FontSize', 22)
        set(gca,'fontsize',22)%, 'XTick', linspace(0,1,6),...
            %'YTick', linspace(0,1,6))
        print('results/Input','-depsc2')
    end

    if ~isempty(triggertime)
        
        sigmaovererr = est_err./err;
        for tau = 1:length(t)
            if any(t(tau)==triggertime)
                sigmaovererr(tau) = scaling_factor;
            end
        end
        
        fig = figure('Units','inches','Position',[0 0 10 4.5],...
            'PaperPositionMode','auto');
        set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
            'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
            'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
            'latex')
        hold on
        plot(t,sigmaovererr, 'b', 'LineWidth', 2);
        plot(triggertime,scaling_factor*ones(size(triggertime)),'ob', 'LineWidth',2)
        %hold on
        %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
        grid off
        axis([0 TEnd -0.1*scaling_factor ...
            2*scaling_factor])
    %     legend({'Adaptive bs.','$$N=5$$',...
    %         '$$N=50$$'},...% '$$N=4500$$'},...
    %         'FontName','CMU Serif');
%         title({'Simulation Results'}, 'FontSize', 22)
        xlabel('Time t (sec)', 'FontSize', 22)
        ylabel('$$\Vert {\sigma}_{\kappa}\Vert / \Vert e \Vert$$', 'FontSize', 22)
        set(gca,'fontsize',22)%, 'XTick', linspace(0,1,6),...
%             'YTick', linspace(0,1,6))
        scale = 0.1;
        pos = get(gca, 'Position');
        pos(2) = pos(2)+scale*pos(4);
        pos(4) = (1-scale)*pos(4);
        set(gca, 'Position', pos)
        MagInset(fig, -1, [triggertime(end-40)*0.98 triggertime(end)*1.02 ...
            0.8*scaling_factor 1.2*scaling_factor], ...
            [triggertime(end)*1.3 9.5 1.2*scaling_factor 1.9*scaling_factor], {'NW','NW';'SE','SE'});        
        ax = gca;
        grid off
        print('results/Trigger','-depsc2')
    end
    
end
