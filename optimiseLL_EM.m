function [Result,param] = optimiseLL_EM(likfun, param, data,nStarts,nsub,parallel,prior)
        
    % initialization
    tol = 1e-3;
    maxiter = 20;
    iter = 0;

    K = length(param);  % Number of parameters
    S = length(data);   % Number of subjects

    if nargin == 7
        m = prior.mean;
        v = prior.var;
    else
        m = randn(1,K);
        v = ones(1,K)*100;   
    end
    
    % identity link function is default
    if ~isfield(param, 'link')
        for k = 1:K
            param(k).link = @(x) x;
        end
    end
    
    
    while iter < maxiter

        iter = iter + 1;
        disp(['.. iteration ',num2str(iter)]);

        % construct prior
        for k = 1:K
            param(k).logpdf = @(x) -0.5 * ((param(k).link(x) - m(k))./sqrt(v(k))).^2 - log((sqrt(2*pi) .* sqrt(v(k))));
        end

        % E-step: Find parameter estimate: 
        if parallel == 1
             Result = optimiseLL_parallel(likfun, param, data,nStarts,nsub);
        else
            Result = optimiseLL(likfun, param, data,nStarts,nsub);
        end
        
        % M-step: update group-level parameters
        v = zeros(1,K);
        for s = 1:S
            v = v + Result.paramfit(s,:).^2 + diag(pinv(Result.H{s}))';
            try
                h = logdet(Result.H{s},'chol');
                L(s) = Result.logpost(s) + 0.5*(Result.K*log(2*pi) - h);
                goodHessian(s) = 1;
            catch
                warning('Hessian is not positive definite');
                try
                    h = logdet(Result.H{s});
                    L(s) = Result.logpost(s) + 0.5*(Result.K*log(2*pi) - h);
                    if isreal(L(s))
                        L(s) = L(s);
                    else
                        L(s) = nan;
                    end
                    goodHessian(s) = 0;
                catch
                    warning('could not calculate');
                    goodHessian(s) = -1;
                    L(s) = nan;
                end
            end
        end
        
        m = nanmean(Result.paramfit);
        v = max(1e-5, v./S - m.^2);  % make sure variances don't get too small
        L(isnan(L)) = nanmean(L); % interpolate to avoid imaginary numbers
        lme(iter) = sum(L) - K*log(sum([data.Nch]));
        
        Result.group.m = m;
        Result.group.v = v;
        Result.lme = lme;
        Result.goodHessian = goodHessian;        
        
        
        if iter > 1 && abs(lme(iter)-lme(iter-1))<tol
            break;
        end

    end
end