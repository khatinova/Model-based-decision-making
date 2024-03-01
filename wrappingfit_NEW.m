
clear
% Set options
optimize = 1; % Set if you wish to optimize free parameters
GeneratePlots = 0; %If you want to generate figures
GenerateSurrData = 0; %Do you want to generate surrogate data?
GenerateCondPlots = 1;
plotDawGraphs = 0;
%% Extracting variables from Gorilla files
filenames= { 
{
'PD101_OFF.mat', 'PD102_OFF.mat', 'PD103_OFF.mat', 'PD104_OFF.mat'...
'PD105_OFF.mat','PD106_OFF.mat','PD107_OFF.mat','PD108_OFF.mat','PD110_OFF.mat'...
'PD111_OFF.mat', 'PD112_OFF.mat', 'PD113_OFF.mat', 'PD114_OFF.mat'...
'PD115_OFF.mat','PD116_OFF.mat','PD117_OFF.mat','PD119_OFF.mat','PD120_OFF.mat'...
}
{
'PD101_ON.mat', 'PD102_ON.mat', 'PD103_ON.mat', 'PD104_ON.mat'...
'PD105_ON.mat','PD106_ON.mat','PD107_ON.mat','PD108_ON.mat','PD110_ON.mat'...
'PD111_ON.mat', 'PD112_ON.mat', 'PD113_ON.mat', 'PD114_ON.mat'...
'PD115_ON.mat','PD116_ON.mat','PD117_ON.mat','PD119_ON.mat','PD120_ON.mat'...
}
{
'C101.mat', 'C102.mat', 'C103.mat', 'C104.mat', 'C105.mat', 'C106.mat', 'C107.mat', ...
'C108.mat', 'C109.mat', 'C110.mat', 'C111.mat', 'C112.mat', 'C113.mat', 'C114.mat',...
'C115.mat', 'C116.mat', 'C117.mat', 'C118.mat' 'C119.mat'
}
};
% initialize the combined group tables for a) trial-wise data (stick models) and b) ratingschange
stayprobTable = table; 

for group = 1:numel(filenames)
    GroupData = []; allratings = []; figure %(for ratinsgchange subplot); %initialize structures and figures 
    for i = 1:length(filenames{group, 1}) % for the number of cells in filenames
        filename = filenames{group}{i}; %extract the filename corresponding to the current subject & group
        result = load(filename); %load the datafile
        %if i == 1; result.data = result.KH290124; end %uncomment if using KHMOCK as the first validation
        if group == 1; Med = 0; Dis = 1; %set PDoff to default group
            if i < 10; subjectData.LearningEffect = 1; % for these participants, this session was the 2nd visit
            else subjectData.LearningEffect = 0; % for these participants, this session was the first
            end
        elseif group == 2; Med = 1; Dis = 1; % set PDon to med = 1, dis = 0 to measure effects of emd
            if i >= 10; subjectData.LearningEffect = 0; % for these participants, this session was the first
            else subjectData.LearningEffect = 1; % for these participants, this session was the 2nd visit
            end
        elseif group == 3; Med = 0; Dis = 0; subjectData.LearningEffect = nan; 
            % set controls to disease = 1 to measure effects of disease
            % no learning effect, as they only did one session
        end
        subjectData.Med = Med;
        subjectData.Dis = Dis;
        % Extract variables for the current subject
        subjectData.group = group;
        subjectData.ID = i; % set subject ID to the current loop number
       
        if group == 3; subjectData.ID = i+18 
        end % for controls, give them new subject IDs
        
        subjectData.A = [result.data.colour_s1, result.data.response_s2]; 
        %action at s1 (red = 1, blue = 2) and shape at s2 (1-4)
        subjectData.Ashape = [result.data.response_s1]; % which s1 shape was chosen (1-5)
        subjectData.R = [result.data.reward]; % did they get a reward
        subjectData.S = [ones(size(result.data.s2_state)), result.data.s2_state]; 
        %state at s1 (always 1) and s2 (state 1 = shapes 1/2, state 2 = shapes 3/4)
        subjectData.shapeIdx = [ floor(result.data.irrels/10), mod(result.data.irrels,10) ]; 
        % what shapes were presented at stage 1 (red | blue)
        
        subjectData.Nch = length(result.data.s3_shape); %number of trials
        subjectData.ratingIdx = [result.data.s3_shape]; %which shape was presented for rating (1-5)?
        subjectData.trans = [result.data.isConsistentMapping]; %common trial (0 = no, 1 = yes)
        subjectData.startTransferChoice = result.data.t_start_s3'; %time rating shape presented
        subjectData.endTransferChoice = result.data.t_response_s3'; %time rating response submitted
        subjectData.start = result.data.t_start_s1'; %s1 time of presentation
        subjectData.finish = result.data.t_response_s1'; %response time on s1
        subjectData.correct = result.data.correct'; % was the s1 choice correct (s2 shape == high value shape)?
        subjectData.LRchoice = result.data.location_s1; % was the chosen s1 colour left or right?
        subjectData.high_value_shape = [result.data.high_value_shape]; %which shape had the high value
 
        %make reward probability matrix for s2 shapes based on high value shape (p(r) = 0.8)
        rewprobValue = 0.2; rewprobArray = zeros(2, 2, 50);
        for k = 1:subjectData.Nch
            highValueIndex = subjectData.high_value_shape(k);
            rewprobArray(:, :, k) = rewprobValue; % Set all values to 0.2
            if highValueIndex == 1; rewprobArray(1, 1, k) = 0.8;
            elseif highValueIndex == 2; rewprobArray(1, 2, k) = 0.8;
            elseif highValueIndex == 3; rewprobArray(2, 1, k) = 0.8;
            else; rewprobArray(2, 2, k) = 0.8;
            end

        end
        subjectData.rewprob = rewprobArray;       
        subjectData.Performance = sum(subjectData.R)/subjectData.Nch; % measure performance

        x = tiedrank(result.data.rating_Coord);% tiedrank the ratings to eliminate subject-specific ranges
        subjectData.zRating = x'; % x - nanmean(x(:,1:10),2); % remove the mean of the 1st 10 trials from the tiedranks
        if result.data.nl_orientation == 1 % what orientation was the number line? (1 = vert, 0 = horizontal)
            subjectData.Rating = (result.data.rating_Coord - 230)'; % subtract the value of coordinate at start of numberline
        else subjectData.Rating = (result.data.rating_Coord - 750)';
        end
        subjectData.Rating(subjectData.Rating < 0) = NaN; 
        % eliminate any ratings that were initiall zero (ie nan), and now will be negative 
    
       % PD Graphing - calculates model-based and model-free behaviour for
       % a) choices 
       % b) ratings
       % c) makes a trialwise table (subjectTable) that has stick probability trial 2:end, based on 
       % consistency, reward, correct, trial + subject ID
        [pMB, pMF, pwin, barv, baraa, MBMFD, subjectTable, pMTbA, pMFbA, pMBTPOSA, pMBNEG, ...
        RTTP, baraaq, ratingschange_chosens] = PDGraphingTR_K(subjectData, i);
        

        % Add the Graphing variables to SubjectData
        subjectData.pMB = pMB;
        subjectData.pMF = pMF;
        subjectData.pMTbA = pMTbA; % model based value transfer (based on last choice of current rated shape win/common)
        subjectData.pMFbA = pMFbA; % model free value transfer (based on last choice of current rated shape win/common)
        subjectData.ratingschange_chosens = ratingschange_chosens;

        stayprobTable = [stayprobTable; subjectTable];
        GroupData = [GroupData;subjectData];
        subjectRatings = subjectData.Rating;
        %allratings = [allratings; subjectRatings];
       
        clear subjectRatings
        %allbarv(i,:)=barv; % barv is the winc/common split of choices
        %allbaraa(i,:)= baraa; %win/common split of shape ratings, chosen not re-seen
    end
    % consolidate/organize the extracted data for each group
    % allratings_g{group} = allratings;
    %GroupTable = struct2table(GroupData);
    %Performance_g{group} = GroupTable.Performance; 
    % insert whichever variable I want to plot a barchart for the three groups at the end of the script

    if plotDawGraphs == 1
      graphTypes = ["Classic Daw", "Ratings Daw"];
        for idx = 1:length(graphTypes)
            figure;   
            if idx == 1
                data = allbarv;
                ylabel('Stay Probability');
                ylim([0.5 0.9]);
            else
                data = allbaraa;
                ylabel('Tiedranked Ratings');
                ylim([0, 250]);
            end
            if group == 1:2
                data(15,:) = nan; % eliminate PD116 participant, as they only chose one colour at s1
            else data(7,:) = nan; % eliminate C107, as they also only chose one colour at s1 (save for one trial)
            end
            errbar = nanstd(data) ./ sqrt(sum(~isnan(data))); %SEM
            x = 1:4;
            barWidth = 0.4;
            hold on;
            bar(x([1, 3]), nanmean(data(:,[1, 3])), barWidth, 'FaceColor', 'red');
            bar(x([2, 4]), nanmean(data(:,[2, 4])), barWidth, 'FaceColor', 'blue');
            er = errorbar(x, nanmean(data), errbar, -errbar);
            er.Color = [0 0 0]; er.LineStyle = 'none';
            xlim([0.5, 4.5]); xticks(x);
            xLabels = ["Win/Common", "Win/Rare", "Loss/Common", "Loss/Rare"];
            xticklabels(xLabels);
            title([graphTypes(idx) ' Graph ' getGroupLabel(group)]);
            hold off
        end
    end
    
    % Fit hybrid learning model
    if optimize == 1
        fit = fitWrapper(GroupData, 'LLmodelRating_K', 1, 0); %data, fucntion, EM, parallel
    else 
        fitall = load('/Users/klarahatinova/Documents/MBDM_2023/Saved_variables/fittedParamCell.mat');
        fit = fitall.ans{group,1}; clear fitall;
    end

    rpeData = [];
    for subject = 1:length(GroupData)
        y = [fit.results.paramfit(subject, :)]; %load the parameters from file for specific subject
        [LL, RPEs] = LLmodelRating_extraction(y, GroupData, subject); % extract RPEs based on fitted parameters
        rpeData(subject).subject = RPEs;
    end

    if GenerateCondPlots == 1
        figure;
        for subject = 1:length(GroupData)
            RC = GroupData(subject).ratingschange_chosens(:);
            MBdQ = rpeData(subject).subject.RPE4(:);
            
            subplot(4,5,subject)
            
            % Plot data
            [Xb,Yb, p,t, h, resid] = conditionalPlot(MBdQ, RC);
            
            % Add dashed lines at y = 0 and x = 0
            hold on;
            line([-1, 1], [0, 0], 'Color', 'k', 'LineStyle', '--');
            line([0, 0], [-50, 50], 'Color', 'k', 'LineStyle', '--');
            hold off;
            
            % Set axis limits
            xlim([-1, 1]);
            ylim([-40, 40]);
            
            title(['MBdQ for Subject ' num2str(subject), ' ', getGroupLabel(group)]);
        end


    end
end

%% Graphs from initial GraphingScriptPDTR
if GeneratePlots == 1;
    figure;
    for i=1:3 % graph of ratings for all three groups
        x = allratings_g{i}(:,2:end); % remove first 2 trials
        % for j=1:size(x,1)
        %     x(j,:) = tiedrank(x(j,:)) / sum(~isnan(x(j,:))) ;
        % end
       % x = x - nanmean(x(:,1:10),2); % subtract off baseline i.e. 1st 10 trials
        errorBarPlot(smoothn(2, x,100),'area',1,'color',colourMap(3,i), 'plotindividuals', 1)
        % remove plotind, 1 to get back to the classic aragement (i,3) 
        hold on
    end
    hold off
    legend({'OFF','','ON','','CTL',''}) %legend
    title 'Ranked ratings over time' % title
    
    
    figure;
    for i = 1:3 % plot whichever PDgraphing variable I hid in Performance_g
        x = nanmean(Performance_g{1,i});
        sem = nanstd(Performance_g{1,i}) / sqrt(sum(~isnan(Performance_g{1,i})));  % Calculate SEM    
        hold on
        bar(i, x, 'FaceColor', colourMap(3, i));
        errorbar(i, x, sem, 'k', 'LineWidth', 1.5);  % Add error bars in black
    end
    hold off
    % Add labels and title
    xlabel('Group'); ylabel('Mean Performance'); legend({'PD off', 'PD on', 'Controls'});
    title('Mean Correct s2 Choices');
    
    figure; %make raw plot of ratings over trials
    for subject = 1:length(GroupData)
        subplot(4, 5, subject)
        plot(GroupData(subject).Rating)
        xlabel('Trial'); ylabel('Raw Rating'); title(['Control #' num2str(subject)])
    end
end
