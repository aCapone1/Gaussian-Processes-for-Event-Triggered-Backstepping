function [] = plot_eventtrig_res(t,x,u,est_err, name, triggertime)
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
    beta = 3;
    scaling_factor = 2*sqrt(beta)*1/Kmin;

    
    TEnd = t(end);
    err1 = abs(x(:,1)-sin(2*pi*t));
    err2 = abs(x(:,2)-2*pi*cos(2*pi*t));
    err3 = abs(x(:,3)-x(:,6));
    errtot = err1+err2+err3;
    
    
%     
%     figure('Units','inches','Position',[0 0 9 4.5],...
%     'PaperPositionMode','auto');
%     set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
%     'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
%     'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
%     'latex')
%     plot(t,err1, 'b', 'LineWidth', 2);
%     %hold on
%     %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
%     grid off
%     axis([0 TEnd -0.1*max(err1) ...
%     1.1*max(err1)])
% %     legend({'Adaptive bs.'},...% '$$N=4500$$'},...
% %     'FontName','CMU Serif');
% %     title({'$$x_1-x_{1,d}$$'}, 'FontSize', 22)
%     xlabel('Time t (sec)', 'FontSize', 22)
%     ylabel('$$e_1$$', 'FontSize', 22)
%     set(gca,'fontsize',22, 'XTick', linspace(0,TEnd,6))
%     print(['Error1',name],'-depsc2')
%     
%     
%     figure('Units','inches','Position',[0 0 9 4.5],...
%     'PaperPositionMode','auto');
%     set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
%     'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
%     'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
%     'latex')
%     plot(t,err2, 'b', 'LineWidth', 2);
%     %hold on
%     %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
%     grid off
%     axis([0 TEnd -0.1*max(err2) ...
%     1.1*max(err2)])
% %     legend({'Adaptive bs.'
% %     'Bs. with GPs, $$N=50$$'},...% '$$N=4500$$'},...
% %     'FontName','CMU Serif');
% %     title({'$$x_2-x_{2,d}$$'}, 'FontSize', 22)
%     xlabel('Time t (sec)', 'FontSize', 22)
%     ylabel('$$e_2$$', 'FontSize', 22)
%     set(gca,'fontsize',22, 'XTick', linspace(0,TEnd,6))
%     print(['Error2',name],'-depsc2')
%     
%     figure('Units','inches','Position',[0 0 9 4.5],...
%     'PaperPositionMode','auto');
%     set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
%     'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
%     'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
%     'latex')
%     plot(t,err3, 'b', 'LineWidth', 2);
%     hold on
%     plot(t,err3, 'r', 'LineWidth', 2);
%     %hold on
%     %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
%     grid off
%     axis([0 TEnd -0.1*max(err3) ...
%     1.1*max(err3)])
% %     legend({'Adaptive bs.'
% %     'Bs. with GPs, $$N=50$$'},...% '$$N=4500$$'},...
% %     'FontName','CMU Serif');
% %     title({'$$x_3-x_{3,d}$$'}, 'FontSize', 22)
%     xlabel('Time t (sec)', 'FontSize', 22)
%     ylabel('$$e_3$$', 'FontSize', 22)
%     set(gca,'fontsize',22, 'XTick', linspace(0,TEnd,6))
%     print(['Error3',name],'-depsc2')

    figure('Units','inches','Position',[0 0 10 4.5],...
    'PaperPositionMode','auto');
    set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
    'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
    'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
    'latex')
    hold on
    plot(t,errtot, 'b', 'LineWidth', 2);
    %hold on
    %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
    grid off
    axis([0 TEnd 0.05 50]) %-0.1*max(errtot) ...
%     1.1*max(errtot)])
%     legend({'Adaptive bs.'
%     'Bs. with GPs, $$N=50$$'},...% '$$N=4500$$'},...
%     'FontName','CMU Serif');
%     title({'$$x-x_d$$'}, 'FontSize', 22)
    xlabel('Time t (sec)', 'FontSize', 22)
    if isequal(name,'triggered')
        ylabel('Tracking err. $$\Vert e\Vert$$', 'FontSize', 22)
    end
    set(gca,'YScale','log')
    set(gca,'fontsize',22)
    set(gca,'XTick', linspace(0,TEnd,6))
    print(['TotalError',name],'-depsc2')
    
    
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
        print(['Input',name],'-depsc2')
    end

    if ~isempty(triggertime)
        
        sigmaovererr = est_err./err;
        for tau = 1:length(t)
            if any(t(tau)==triggertime)
                sigmaovererr(tau) = max(sigmaovererr);
            end
        end
        
        est_err = est_err/max(est_err)*0.1;
        fig = figure('Units','inches','Position',[0 0 9 4.5],...
            'PaperPositionMode','auto');
        set(0,'defaulttextinterpreter','latex','DefaultTextFontname', ...
            'CMU Serif','DefaultTextFontSize',22,'DefaultLegendFontname', ...
            'CMU Serif','DefaultLegendFontSize',22,'defaultlegendinterpreter',...
            'latex')
        hold on
        plot(t,sigmaovererr, 'r', 'LineWidth', 2);
        plot(triggertime,max(sigmaovererr)*ones(size(triggertime)),'or', 'LineWidth',2)
        %hold on
        %plot(t3, x3(:,1), 'c--', 'LineWidth', 2);
        grid off
        axis([0 TEnd -0.1*max(sigmaovererr) ...
            1.5*max(sigmaovererr)])
    %     legend({'Adaptive bs.','$$N=5$$',...
    %         '$$N=50$$'},...% '$$N=4500$$'},...
    %         'FontName','CMU Serif');
%         title({'Simulation Results'}, 'FontSize', 22)
        xlabel('Time t (sec)', 'FontSize', 22)
        ylabel('$$\Vert \bar{\sigma}\Vert / \Vert e \Vert$$', 'FontSize', 22)
        set(gca,'fontsize',22)%, 'XTick', linspace(0,1,6),...
%             'YTick', linspace(0,1,6))
        MagInset(fig, -1, [2.5 3.2 700 900], [4.5 9.5 200 1200], {'NW','NW';'SW','SW'});        ax = gca;
        grid off
        print(['Trigger',name],'-depsc2')
    end
    
end
