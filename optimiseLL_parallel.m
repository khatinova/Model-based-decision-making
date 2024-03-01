function Result = optimiseLL_parallel(likfun, param, data,nStarts,nsub)

    options = optimset('Display','off','MaxFunEvals',2000);
    warning off all
        
    K = length(param);  % Number of parameters
    
    Result.K = K;
    Result.S = nsub;

    % save info to results structure
    Result.param = param;
    Result.likfun = likfun;
    
    % preallocate
    [logpost, logL, BIC, AIC, exFlags] = deal(nan(nsub, 1));
    [x_out] = nan(nsub, K);
    [H] = cell(nsub, 1);

    parfor s = 1:nsub
        
        disp(['Subject ',num2str(s)]);
        
        % construct posterior function
        fpost = @(x) -likfun_post(x,param,data(s),likfun);
        
        for i = 1:nStarts
            success = 0;
            x0 =  randn(1,K); %[[1.5 0.9 0.1 -0.2 1 1 0.2]'+.2*randn(7,1)]'
            while ~success
                try
                    [x,nlogp,exFlag,~,~,h] = fminunc(fpost,x0,options);
                    success = 1;
                catch 
                    x0 = randn(1,K);
                end
            end
            
            logp = -nlogp;
            if i == 1 || logpost(s) < logp
                logpost(s) = logp;
                logL(s) = likfun(x,data(s));
                x_out(s,:) = x;
                H{s} = h;
                exFlags(s) = exFlag;
            end
            if i > 20
                if exFlags(s) == 1 || exFlags(s) == 2 || exFlags(s) == 3
                    break
                end
            end
        end
        
       BIC(s) = K*log(data(s).Nch) - 2*logL(s);
       AIC(s) = K*2 - 2*logL(s);
        
    end
      
    % save
    Result.logpost = logpost';
    Result.logL = logL;
    Result.BIC = BIC;
    Result.AIC = AIC;
    Result.paramfit = x_out;
    Result.H = H;
    Result.exFlag = exFlags;

end


