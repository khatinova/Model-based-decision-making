function logp = likfun_post(x,param,data,likfun)
    
    % Evaluate log probability of parameters under the (unnormalized) posterior.
    %
    % USAGE: logp = likfun_post(x,param,data,likfun)
    %
    % INPUTS:
    %   x - parameter values
    %   param - parameter structure
    %   data - data structure
    %   likfun - function handle for likelihood function
    %
    % OUTPUTS:
    %   logp - log unnormalized posterior probability
    %

    logp = likfun(x,data);
    
%     xtrans = [exp(x(1)) exp(x(2)) exp(x(3)) 1./(1+exp(x(4))) 1./(1+exp(x(5))) 1./(1+exp(x(6))) x(7)];
    if isfield(param,'logpdf')
        for k = 1:length(param)
            logp = logp + param(k).logpdf(x(:,k));
        end
    end
end