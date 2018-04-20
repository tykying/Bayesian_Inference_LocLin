function ylim_vis = obtain_ylim_OrdMag(MaxData, MinData)
        base = 10;
        
        ylim_ub = max([abs(MaxData), abs(MinData)]);
        
        base = 10;
        % for a real number = L*10*M, evaluate L and M. L between 0 and 10
        Ub_magorder = floor(log(ylim_ub)./log(base));
        Ub_lead = ylim_ub/(10^Ub_magorder);
        
        if Ub_lead >= 5
            ylim_vis_b = base^(Ub_magorder+1);
%         elseif Ub_lead < 2.5
%             ylim_vis = 0.25*10^(Ub_magorder+1);
        else
            ylim_vis_b = 0.5*base^(Ub_magorder+1);
        end
        
        % Round to nearest 5
        if sign(MaxData) * sign(MinData) < 0
            % Opposite sign
            ylim_vis = [-ylim_vis_b, ylim_vis_b];
        else
            % Opposite sign
            if sign(MaxData) > 0
                ylim_vis = [0, ylim_vis_b];
            else
                ylim_vis = [-ylim_vis_b, 0];
            end
        end
end