function [l2err_triggered, l2err_untriggered, Ndata_triggered] =  compute_l2err(nreps)

    l2err_triggered = [];
    l2err_untriggered = [];
    Ndata_triggered = [];
    
    for rep =1:nreps
        
        load("results/results_evnttrig" + int2str(rep))
        try
            time_untrig = t_untrig;
            state_untrig = x_untrig;
        catch
        end
        
        %Desired signals and derivatives, triggered
        x1d = sin(2*pi*time);
        z11 = state(:,4);
        z21 = state(:,6);
        x2d = z11;
        x3d = z21;

        err2_trig = ...
            ((state(:,1)-x1d).^2 + (state(:,2)-x2d).^2 + (state(:,3)-x3d).^2);
        l2err_triggered(rep) = sqrt(trapz(time,err2_trig));
        Ndata_triggered(rep) = length(Xd1);
        
        %Desired signals and derivatives, untriggered
        x1d = sin(2*pi*time_untrig);
        z11 = state_untrig(:,4);
        z21 = state_untrig(:,6);
        x2d = z11;
        x3d = z21;

        err2_untrig = ...
            ((state_untrig(:,1)-x1d).^2 + ...cle
            (state_untrig(:,2)-x2d).^2 + (state_untrig(:,3)-x3d).^2);
        l2err_untriggered(rep) = sqrt(trapz(time_untrig,err2_untrig));
    end

end