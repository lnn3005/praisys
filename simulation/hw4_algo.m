

C = 30;         % fixed cost of repair
bbeta = 0.5;    % rate of deterioration
ggamma = 0.6;   % discount factor
xmax = 10;      %





% construct the value function approximation
pphi = {@(x) x.^0, @(x) x.^1, @(x) x.^2, @(x) x.^3};
n_basis = numel(pphi);
ttheta = zeros(n_basis,1);

itermax = 20;
xx_plot = linspace(0,xmax,101);
show_samples = true;
cvx_solver SeDuMi
cvx_quiet true
for iter = 1:itermax
    
    % sample N basepoint states uniformly in [0, xmax]
    N = 30;
    xx = rand(N,1)*xmax;
    
    % for each base state, sample next state and reward
    M = 30;
    xxnext = zeros(N,M,2);
    rewards = zeros(N,M,2);
    VVxxnext = zeros(N,M,2);
    for ia = 1:2
        a = ia-1;
        for ii = 1:N
            x = xx(ii);
            for jj = 1:M
                [xxnext(ii,jj,ia) rewards(ii,jj,ia)] = ...
                    myReplacementProblemEnvironment(x,a, bbeta, C, xmax);
            end
        end
    end
    PPhi_sum = zeros(N,M,2);
    for ii = 1:n_basis
        PPhi_ii = ttheta(ii)*pphi{ii}(xxnext);
        PPhi_sum = PPhi_sum+PPhi_ii;
    end
    QQ_sample = rewards + ggamma*PPhi_sum;
    QQ_mean = squeeze(sum(QQ_sample,2)/M);
    Vhat_xx = max(QQ_mean,[],2);
    
    % update the value function approximation by solving a regression problem
    PPhi_x = zeros(N,n_basis);
    for ii = 1:n_basis
        PPhi_x(:,ii) = pphi{ii}(xx);
    end
    
    cvx_begin
        variables ttheta_new(n_basis,1)
        variables Vmodel_xx(N,1);
        minimize(norm(Vmodel_xx-Vhat_xx,2))

        subject to
        Vmodel_xx == PPhi_x*ttheta_new

        %don't penalize the theta of the constant term
        norm(ttheta_new(2:end),2) <= 20; 
    
    cvx_end
    
    
    
    ttheta = ttheta_new;
    
    % plot the new value function
    
    %{
    PPhi_plot = zeros(numel(xx_plot),n_basis);
    for ii = 1:n_basis
        PPhi_plot(:,ii) = pphi{ii}(xx_plot);
    end
    VV_plot = PPhi_plot*ttheta_new;
    if iter == 1
        figure;
        hold on
        if show_samples
            hdl_samples = plot(xx,Vhat_xx,'.');
        end
        hdl = plot(xx_plot,VV_plot);
        hdl_t = title('iter = 1');
        xlabel('state x')
        ylabel('value function approximation')
    else
        set(hdl,'ydata',VV_plot);
        if show_samples
            set(hdl_samples,'xdata',xx,'ydata',Vhat_xx);
        end
        set(hdl_t,'string',['iter = ' num2str(iter)]);
        drawnow
    end
    %}
end