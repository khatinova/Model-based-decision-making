function Result = optimiseLL(likfun, param, data,nStarts,nsub)

options = optimoptions('fminunc', 'Display', 'off', 'MaxFunEvals', 2000, 'SpecifyObjectiveGradient', true);

warning off all

    for sub = 1:nsub
        
        disp(['Subject ',num2str(sub)]);

        % construct posterior function
        fpost = @(x) -likfun_post(x,param,data(sub),likfun);

        K = length(param);  % Number of parameters

        for n = 1:nStarts
            success = 0;
            x0 =  randn(1,K); %[[1.5 0.9 0.1 -0.2 1 1 0.2]'+.2*randn(7,1)]'
%             while ~success
%                 try
                    [x,nlogp,exFlag,~,~,H] = fminunc(fpost,x0,options);
%                     success = 1;
%                 catch 
%                     x0 = randn(1,K);
%                 end
%             end
            
            logp = -nlogp;
            if n == 1 || Result.logpost(sub,1) < logp
                Result.logpost(sub,1) = logp;
                Result.LogL(sub,1) = likfun(x,data(sub));
                [~, Result.N(sub,1)] = likfun(x,data(sub));
                Result.paramfit(sub,:) = x;
                Result.H{sub} = H;
                Result.exFlag(sub,1) = exFlag;
            end
            Result.BIC(sub,1) = K*log(Result.N(sub)) - 2*Result.LogL(sub,1); % 402 trials in total (201 per session)
            Result.AIC(sub,1) = K*2 - 2*Result.LogL(sub,1);
            Result.param = param;
            Result.K = K;
            if n > 20
                if Result.exFlag(sub) == 1 || Result.exFlag(sub) == 2 || Result.exFlag(sub) == 3
                    break
                end
            end
        end
    end
end


