function [xnext, reward] = myReplacementProblemEnvironment(x,a,bbeta,C, xmax)
% bbeta: rate of deterioration. C: repair cost.
% the maintenance cost is c(x)=log(1+x).
% convention: a==0 = keep, a == 1 = replace
% automatic replacement if y > xmax.
    % draw (y-x) following exponential of rate bbeta
    
  
    
    if a == 0 %keep
        xnext = x-log(rand)/bbeta; % deterioration is exponentially distributed
        if xnext > xmax
            % replace automatically
            xnext = -log(rand)/bbeta;
            c0 = 4*0;
            reward = -C-c0;
        else
            cx = 4*x;
            reward = -cx;
        end
    else % replace
        xnext = -log(rand)/bbeta;
        c0 = 4*0;
        reward = -C-c0;    
    end