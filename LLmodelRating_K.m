function [LL, N, dsurr] = LLmodelRating_K(y,D,opts) 

% Find parameters estimates for original reduced model with 5
% parameters. Learning rates and softmax parameteres for stage 1 and 2 
% are assumed to be the same.
% Transition probabilities are associated with colour. Shapes presented
% at stage 1 are irrelevant, but have acquired value measured through 
% rating. Optimise model to account for this.
%
% USAGE: [LL, N,dsurr] = LLmodelRating(y,data,opts)
%
% INPUTS:
% Free parameters (y):
% b1          = y(1);
% b2          = y(1);
% alpha 1     = y(2);
% alpha 2     = y(2);
% lambda      = y(3);
% omega       = y(4);
% pi          = y(5);
% data = subject line data
% opts = generate surrogate data?

%
% OUTPUTS:
%   LL          = Log likelihood 
%   N           = number of trials
%   dsurr       = latent variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = 192;
    LL = 0;
    
    b1      = exp(y(1));
    b2      = b1;
    alpha1  = 1./(1+exp(-y(2)));
    alpha2  = alpha1;
    lambda  = 1./(1+exp(-y(3)));
    w       = 1./(1+exp(-y(4)));
    st      = exp(y(5));

    if opts.generatesurrogatedata == 1
        rewprob = D.rewprob; %reward probability matrix
        trans = D.trans; %common or rare transition
    end   

        for t = 1:N

            s1 = 1;                      % stage 1 state (A)
            s2 = D.S(t, 2) + 1;          % stage 2 state (A) % which state you got to
            c1 = D.A(t, 1);              % colour choice stage 1
            c2 = D.A(t, 2);              % shape at s2
            r  = D.R(t);                 % reward
            if s2 == 2 
                c2 = c2 - 2;
            end
             % choice stage 2 (must be 1 or 2) CHANGED FROM ORIGINAL

            if t == 1 
                Qd = zeros(3,2);            % Q(s,a): state-action value function for Q-learning
                M = [0 0];
                count = zeros(2);
                Qr = [0 0];
            end   
            
            % Break if trial was missed
            if (isnan(c1) || isnan(c2)) && opts.generatesurrogatedata == 0
                continue
            end % in my code, I've used Nan to denote missed trials
        

            if count(1,1) + count(2,2) > count(1,2) + count(2,1)        % find most likely transition
			    Tr = .3+.4*eye(2);                                    
                % creates matrix where diagonal values are 0.7 and others are 0.3
            else
			    Tr = .3+.4*(1-eye(2));
                % creates matrix where diagonal values are 0.3 and others are 0.7
            end
            
            maxQ = max(Qd(2:3,:),[],2);                            % optimal reward at second step
            % determined by extracting columns 2 and 3 from Qd (state action value)
            % [], 2 implies action should be carried out columnwise (according to 2nd dimension
            
            Qm = (Tr'*maxQ)'; %model-based value is the max value combined with transition
    
            Q = w * Qm + (1-w) * Qd(s1,:) + st * M;               % mix TD and model-based values
    
            % first stage choice
            if opts.generatesurrogatedata == 1                           % if surrogate data is generated, generate choice 1
                c1 = (rand > exp(b1*Q(1))/sum(exp(b1*Q))) +1;     
                % make choice using softmax
                if trans(t)
                    if c1 == 1; s2 = 2; % find stage 2 state based on transition
                    else; s2 = 3; 
                    end   
                    else; if c1 == 2; s2 = 3; else; s2 = 2; end
                end
            else
                LL = LL + b1*Q(c1) - log(sum(exp(b1*Q)))                  % update likelihoods stage 1
            end
            
            % second stage choice
            if opts.generatesurrogatedata == 1                              % if surrogate data is generated, generate choice 2
                c2 =  (rand > exp(b2*Qd(s2,1))/sum(exp(b2*Qd(s2,:))))  +1;  % make choice using softmax 
                r = rand < rewprob(s2-1,c2,t);
            else
                LL = LL + b2*Qd(s2,c2) - log(sum(exp(b2*Qd(s2,:))))        % update likelihoods stage 2
            end
            
            dtQ = [];
            dtQ(1) = Qd(s2,c2) - Qd(s1,c1);                         % backup with actual choice (i.e., sarsa)
            Qd(s1,c1) = Qd(s1,c1) + alpha1 * dtQ(1);                % update TD value function
    
            dtQ(2) = r - Qd(s2,c2);                                 % prediction error (2nd choice)
    
            Qd(s2,c2) = Qd(s2,c2) + alpha2 * dtQ(2);                % update TD value function
            Qd(s1,c1) = Qd(s1,c1) + lambda * alpha1 * dtQ(2);       % eligibility trace
    
            if trans(t) == 0 % common
                RPE11 =(r - Qr(c1)); 
                Qr(c1) = Qr(c1) + alpha1 * (r - Qr(c1));
            else % rare
                RPE11 = (r - Qr(3 - c1));
                Qr(3-c1) = Qr(3-c1) + alpha1 * (r - Qr(3- c1));
            end
    
            count(s2-1,c1) =  count(s2-1,c1)+1;
            
            M = [0 0];
            M(c1) = 1;                                              % make the last choice sticky
    
            if opts.generatesurrogatedata == 1                      % store latents
                MB_new = (Tr'* (max(Qd(2:3,:),[],2)))';             
                if mean(rewprob(s2-1,:,t)) > 0.2                     % optimal reward at second step
                    winState = 2;
                else
                    winState = 3;
                end
                dsurr.A(t,1) = c1; dsurr.A(t,2) = c2;               % store choices
                dsurr.S(t,1) = 1; dsurr.S(t,2) = s2-1;              % store states
                dsurr.R(t,1) = r;                                   % store rewards
                dsurr.trans(t,1) = trans(t);                        % store transitions
                dsurr.Qmb_chosen(t,1) = Qm(c1);                     % MB value Stage 1 chosen option       
                dsurr.Qmb_unchosen(t,1) = Qm(3-c1);                 % MB value Stage 1 unchosen option
                dsurr.Qmf_chosen(t,1) = Qd(1,c1);                   % MF value Stage 1 chosen option
                dsurr.Qmf_unchosen(t,1) = Qd(1,3-c1);               % MB value Stage 1 unchosen option
                dsurr.Qmb_new_chosen(t,1) = MB_new(c1);             % new MB value Stage 1 chosen option after update stage 2 (value used on t+1)
                dsurr.Qmb_new_unchosen(t,1) = MB_new(3-c1);         % new MB value Stage 1 unchosen option after update stage 2 (value used on t+1)
                dsurr.Qs2_chosen(t,1) = Qd(s2,c2);                  % value of stage 2 choice
                dsurr.Qs2_unchosen(t,1) = Qd(s2,3-c2);              % value of stage 2 option not chosen
                dsurr.winState(t,1) = winState;                     % state with highest prob to obtain reward
                dsurr.bestStateChosen(t,1) = (winState == s2);      % did participant end up in best state
                if trans(t) == 1 && (winState == s2)
                    dsurr.bestS1choice(t,1) = c1;                   % if transition was common and participant ended up in winstate, best choice is c1.
                else
                    dsurr.bestS1choice(t,1) = 3-c1;                 % else, choice is the other option
                end            
                dsurr.RPE1(t,1) = r - Qm(c1);                       % RPE MB
                dsurr.RPE2(t,1) = r - Qd(1,c1);                     % RPE MF stage 1
                dsurr.RPE3(t,1) = Qd(s2,c2) - Qd(1,c1);             % RPE MF (model-based) - used in original model
                dsurr.RPE4(t,1) = r - MB_new(c1);                   % RPE MB "retrospective"
                dsurr.RPE5(t,1) = r - Qd(s2,c2);                    % RPE MF stage 2 - used in original model
                dsurr.RPE6(t,1) = (Qd(s2,c2) - Qd(1,c1)) +  (r - Qd(1,c1)); % RPE MF combined
                dsurr.RPE7(t,1) = r - Qm(dsurr.bestS1choice(t,1));  % RPE MB that would have led to win state
                dsurr.RPE8(t,1) = r - Qd(1,dsurr.bestS1choice(t,1));% RPE MF that would have led to win state
                dsurr.RPE9(t,1) = r - Qm(s2-1);                     % RPE MB of choice that had greatest likelihood of leading to current state
                dsurr.RPE10(t,1) = r - Qd(1,s2-1);                  % RPE MB of choice that had greatest likelihood of leading to current state
                dsurr.RPE11(t,1) = RPE11;                           % RPE10*
                dsurr.TruePar(1,:) = y;                             % Parameters used 
                dsurr.Nch(1) = N;                                   % number of trials
            end
        
        end
    
end

