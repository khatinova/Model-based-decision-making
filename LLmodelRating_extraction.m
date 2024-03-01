function [LL, RPEs] = LLmodelRating_extraction(y,D, subject) %opts has been reduced to [1 or 0]

% Find parameters estimates for original reduced model with 5
% parameters. Learning rates and softmax parameteres for stage 1 and 2 
% are assumed to be the same.
% Transition probabilities are associated with colour. Shapes presented
% at stage 1 are irrelevant, but have acquired value measured through 
% rating. Optimise model to account for this.
%
% USAGE: [LL, N,RPEs] = LLmodelRating(y,data,opts)
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

%
% OUTPUTS:
%   LL          = Log likelihood 
%   N           = number of trials
%   RPEs       = latent variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = D(subject).Nch;
    LL = 0;
    
    b1      = exp(y(1));
    b2      = b1;
    alpha1  = 1./(1+exp(-y(2)));
    alpha2  = alpha1;
    lambda  = 1./(1+exp(-y(3)));
    w       = 1./(1+exp(-y(4)));
    st      = exp(y(5));


         rewprob = D(subject).rewprob; %reward probability matrix
         trans = D(subject).trans; %common or rare transition


        for t = 1:N

            s1 = D(subject).S(t, 1);              % stage 1 state (always 1)
            s2 = D(subject).S(t, 2) + 1;          % stage 2 state (2 or 3)
            c1 = D(subject).A(t, 1);              % colour choice stage 1 (1 or 2)
            c2 = D(subject).A(t, 2);              % choice stage 2, 1 or 2 (shape determined by s2, c2 combo)
            
            if s2 == 3
                c2 = c2 - 2;
            end
             % choice stage 2 (must be 1 or 2) CHANGED FROM ORIGINAL

            r  = D(subject).R(t);           % reward (0 or 1)
    
            if t == 1 
                Qd = zeros(3,2);            % Q(s,a): state-action value function for Q-learning
                M = [0 0];                  % Identifies which c1 was chosen on previous trial and makes it "sticky"
                count = zeros(2);           % counts which choice was made on previous trial
                Qr = [0 0];                 % Transition-based values (??)
            end   
            
            % Break if trial was missed
            if (isnan(c1) || isnan(c2)) 
                continue
            end % I've used Nan to denote missed trials
        

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
            
            Qm = (Tr'*maxQ)'; %model-based value

            % mix TD and model-based values, st is the choice stickiness, which builds on the previous c1 (M)    
            Q = w * Qm + (1-w) * Qd(s1,:) + st * M;                                                  
    
            % LL of first stage choice
            LL = LL + b1*Q(c1) - log(sum(exp(b1*Q)));             % update likelihoods stage 1
            
            % LL of second stage choice
            LL = LL + b2*Qd(s2,c2) - log(sum(exp(b2*Qd(s2,:))));  % update likelihoods stage 2


            dtQ = [];                                              % temporal difference variable
            dtQ(1) = Qd(s2,c2) - Qd(s1,c1);                        % backup with actual choice (i.e. sarsa)
            Qd(s1,c1) = Qd(s1,c1) + alpha1 * dtQ(1);                % update TD value function
    
            dtQ(2) = r - Qd(s2,c2);                                 % prediction error (2nd choice)
    
            Qd(s2,c2) = Qd(s2,c2) + alpha2 * dtQ(2);                % update TD value function
            Qd(s1,c1) = Qd(s1,c1) + lambda * alpha1 * dtQ(2);       % eligibility trace
    
            if trans(t) == 1 % common
                RPE11 =(r - Qr(c1)); 
                Qr(c1) = Qr(c1) + alpha1 * (r - Qr(c1)); % updates c1 value based on common or rare transitions
            else % rare
                RPE11 = (r - Qr(3 - c1));
                Qr(3-c1) = Qr(3-c1) + alpha1 * (r - Qr(3- c1));
            end
    
            count(s2-1,c1) =  count(s2-1,c1)+1;
            
            M = [0 0];
            M(c1) = 1;               % make the last choice sticky
    
         
            MB_new = (Tr'* (max(Qd(2:3,:),[],2)))';             
            if mean(rewprob(1,:,t)) > 0.2           % optimal reward STATE at second step
                winState = 2;
            else
                winState = 3;
            end
            RPEs.Qmb_chosen(t,1) = Qm(c1);                     % MB value Stage 1 chosen option       
            RPEs.Qmb_unchosen(t,1) = Qm(3-c1);                 % MB value Stage 1 unchosen option
            RPEs.Qmf_chosen(t,1) = Qd(1,c1);                   % MF value Stage 1 chosen option
            RPEs.Qmf_unchosen(t,1) = Qd(1,3-c1);               % MB value Stage 1 unchosen option
            RPEs.Qmb_new_chosen(t,1) = MB_new(c1);             % new MB value Stage 1 chosen option after update stage 2 (value used on t+1)
            RPEs.Qmb_new_unchosen(t,1) = MB_new(3-c1);         % new MB value Stage 1 unchosen option after update stage 2 (value used on t+1)
            RPEs.Qs2_chosen(t,1) = Qd(s2,c2);                  % value of stage 2 choice
            RPEs.Qs2_unchosen(t,1) = Qd(s2,3-c2);              % value of stage 2 option not chosen
            RPEs.winState(t,1) = winState;                     % state with highest prob to obtain reward
            RPEs.bestStateChosen(t,1) = (winState == s2);      % did participant end up in best state
            if trans(t) == 1 && (winState == s2) || (trans(t) == 0 && winState ~= s2)
                RPEs.bestS1choice(t,1) = c1;                   % if transition was common and participant ended up in winstate, best choice is c1.
            else
                RPEs.bestS1choice(t,1) = 3-c1;                 % else, choice is the other option
            end  


            RPEs.RPE1(t,1) = r - Qm(c1);                       % RPE MB
            RPEs.RPE2(t,1) = r - Qd(1,c1);                     % RPE MF stage 1
            RPEs.RPE3(t,1) = Qd(s2,c2) - Qd(1,c1);             % RPE MF (model-based) - used in original model
            RPEs.RPE4(t,1) = r - MB_new(c1);                   % RPE MB "retrospective" ie the value of the chosen s1 colour computed based on max s2 value
            RPEs.RPE5(t,1) = r - Qd(s2,c2);                    % RPE MF stage 2 - used in original model
            RPEs.RPE6(t,1) = (Qd(s2,c2) - Qd(1,c1)) +  (r - Qd(1,c1)); % RPE MF combined
            RPEs.RPE7(t,1) = r - Qm(RPEs.bestS1choice(t,1));  % RPE MB that would have led to win state
            RPEs.RPE8(t,1) = r - Qd(1,RPEs.bestS1choice(t,1));% RPE MF that would have led to win state
            RPEs.RPE9(t,1) = r - Qm(s2-1);                     % RPE MB of choice that had greatest likelihood of leading to current state
            RPEs.RPE10(t,1) = r - Qd(1,s2-1);                  % RPE MB of choice that had greatest likelihood of leading to current state
            RPEs.RPE11(t,1) = RPE11;                           % RPE10*
            RPEs.Par(1,:) = y;                                 % Parameters used 
            RPEs.Nch(1) = N;                                   % number of trials
            
        
        end
    
end

