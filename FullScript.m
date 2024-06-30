%% This is the compiled script designed to:
% 1. extract task and demographic data,
% 2. analyse behaviour and ratings on the task using PDGraphingTR_K
% 3. make a results table for linear models & generate tables to assess the
% significant effects/interactions and model fits
% 4. optimize learning parameters using fitWrapper and LLmodelRating_K (EM works)
% 5. generate RPEs based on the modelled s1 colour / s2 shape value
% 6. look at the best fitting RPE to explain ratings change

%% Set options
clear
options.savetables = 0; % do you want to save tables for the group data, ratings, results tables (short and trialwise)
options.GenerateSurrData = 0; % generate surrogate data?
options.fitRegression = 0; % do you want to regress the ratings against the generated RPEs?
options.GenerateCondPlots = 0; %to be added, but do you want to generate RPE vs ratingschange plots? 
options.optimize = 0; % do you want to fit a RL model to the data using LLmodelRating_K?

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
'C115.mat', 'C116.mat', 'C117.mat', 'C118.mat' 'C119.mat'...
}
};

% initialize the structures and tables 
stayprobTable = table; HCstayprobTable = table; HCGroupData = []; shortTable = []; clear subjectData
groupRPEt = table(); allRPEt = table(); GroupTable = table(); params = [];

PDDemographicData = load("PDDemographicData.mat"); 
HCDemographicData = load("HCDemographicData.mat"); 
f = figure('Name', 'Selected subplots', 'PaperOrientation','landscape');
for group = 1:3
    GroupData = []; groupRPEt = table(); GroupTable = table(); allratings = table();%initialize structures and figures 
    for i = 1:length(filenames{group, 1}) % for the number of cells in filenames
        filename = filenames{group}{i}; %extract the filename corresponding to the current subject & group
        result = load(filename); %load the datafile
        subjectData.ID = i; % set subject ID to the current loop number
        
        if group < 3
            dd = PDDemographicData.data(i,:);
            subjectData.AMI = dd.AMITotal;
            subjectData.AMIB = dd.AMIBehaviour;
            subjectData.AMIS = dd.AMISocial;
            subjectData.AMIE = dd.AMIEmotional;
            subjectData.HADS = dd.HADSTotal;
            subjectData.UPDRST = dd.UPDRSTotal;
            subjectData.Grit = dd.GritTotal; 
            subjectData.Impulsivity = dd.SUPPSP;
            subjectData.Age = dd.Age;
            subjectData.Sex = dd.Sex;
        elseif group == 3
            dd = HCDemographicData.data(i,:);
            subjectData.AMI = dd.AMITotal;
            subjectData.AMIB = dd.AMIBehaviour;
            subjectData.AMIS = dd.AMISocial;
            subjectData.AMIE = dd.AMIEmotional;
            subjectData.HADS = dd.HADSTotal;
            subjectData.UPDRST = 0;
            subjectData.Grit = dd.GritTotal; 
            subjectData.Impulsivity = dd.SUPPSP;
            subjectData.Age = dd.Age;
            subjectData.Sex = dd.Sex;
        end
        if group == 1
            Med = -1; Dis = 1; 
            if i < 10; subjectData.LearningEffect = 1; % for these participants, this session was the 2nd visit
            else; subjectData.LearningEffect = -1; % for these participants, this session was the first
            end
        elseif group == 2
            Med = 1; Dis = 1; % set PDon to med = 1, dis = 1 to measure effects of medication
            if i >= 10 
                subjectData.LearningEffect = -1; % for these participants, this session was the first
            else 
                subjectData.LearningEffect = 1; % for these participants, this session was the 2nd visit
            end
        elseif group == 3
            Med = 0; Dis = 0; %set PDoff to default group = 0 and med = 0
            subjectData.LearningEffect = -1; % no learning effect, as they only did one session
            subjectData.ID = i+18; 
        end

        subjectData.Med = Med;
        subjectData.Dis = Dis;
        % Extract variables for the current subject
        subjectData.group = group;
        subjectData.ID = i; % set subject ID to the current loop number

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
        subjectData.startTransferChoice = [result.data.t_start_s3]; %time rating shape presented
        subjectData.endTransferChoice = [result.data.t_response_s3]; %time rating response submitted
        subjectData.start = [result.data.t_start_s1]; %s1 time of presentation
        subjectData.finish = [result.data.t_response_s1]; %response time on s1
        subjectData.correct = [result.data.correct]; % was the s1 choice correct (s2 shape == high value shape)?
        subjectData.LRchoice = [result.data.location_s1]; % was the chosen s1 colour left or right?
        subjectData.high_value_shape = [result.data.high_value_shape]; %which shape had the high value
 
        subjectData.RT = subjectData.finish - subjectData.start;
        subjectData.version = [result.data.ExperimentVersion];

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
        subjectData.ScreenSize = result.data.ParticipantViewportSize;

        subjectData.nl_orientation = result.data.nl_orientation; nl = subjectData.nl_orientation(1); % what orientation was the number line? (1 = vert, 0 = horizontal)
        
        % Handle repeated rating coordinates - these were trials on which the participant did not respond 
        % and were thus given exactly the same coordinate as the previous trial 
        rating_coord = result.data.rating_Coord;
        nan_idx = false(size(subjectData.RT));
        for t = 2:length(rating_coord)
            if subjectData.RT(t) >= 4000
                nan_idx(t) = true;
            end
        end
        
        % Apply NaN indices to all fields of the missed trial
        subjectData.A(nan_idx, :) = NaN;
        subjectData.Ashape(nan_idx) = NaN;
        subjectData.R(nan_idx) = NaN;
        subjectData.S(nan_idx, :) = NaN;
        subjectData.shapeIdx(nan_idx, :) = NaN;
        subjectData.ratingIdx(nan_idx) = NaN;
        subjectData.trans(nan_idx) = NaN;
        subjectData.startTransferChoice(nan_idx) = NaN;
        subjectData.endTransferChoice(nan_idx) = NaN;
        subjectData.start(nan_idx) = NaN;
        subjectData.finish(nan_idx) = NaN;
        subjectData.correct(nan_idx) = NaN;
        subjectData.LRchoice(nan_idx) = NaN;
        subjectData.high_value_shape(nan_idx) = NaN;
        
        % Initialize the adjusted_rating_Coord array
        adjusted_rating_Coord = nan(size(rating_coord));
        
        % Iterate over each row
        for j = 1:length(rating_coord)
            % Get the screen height and width for the current trial
            xysize = strsplit(string(subjectData.ScreenSize(j)), "x"); % splits the screensize into height and width
            width = str2double(xysize(1));
            height = str2double(xysize(2));
            
            if subjectData.nl_orientation(j) == 1 % if vertical
                % Invert the Y coordinate and normalize by height
                adjusted_rating_Coord(j) = (height - rating_coord(j)) / height;
            elseif subjectData.nl_orientation(j) == 0 % if horizontal
                % Normalize the X coordinate by width
                adjusted_rating_Coord(j) = rating_coord(j) / width;
            end
        end
        
        % Assign the adjusted ratings back to the subjectData structure
       subjectData.adj_coord = adjusted_rating_Coord;
       subjectData.Rating = adjusted_rating_Coord;
       subjectData.Rating(nan_idx | (0.2 > adjusted_rating_Coord > 0.8)) = NaN;
       subjectData.zRating = nanzscore(subjectData.Rating);
       kk(i, group) = nanmean(subjectData.Rating(1:10));
       subjectData.kRating = subjectData.Rating - kk(i,group);


       
       
       %subjectData.zRating = nanzscore(rating_coord);

       [var, subjectTable, ratingschange_chosens] = PDGraphingTR_K1(subjectData, group);
       % Includes calculation of previous previous time seen + difference in ratings 1 trial back
       
       %[var, subjectTable, ratingschange_chosens] = PDGraphingTR_K(subjectData, group); % alternatively 
       
       % PD Graphing - calculates model-based and model-free behaviour for
       % a) choices 
       % b) ratings
       % c) makes a trialwise table (subjectTable) that has stick probability trial 2:end, based on 
       % consistency, reward, correct, trial + subject ID
       % Has two versions: PDGraphing_K is the longform, PDGraphing_K1 is shortform (maybe more accurate)
   
        % Add the Graphing variables to subjectData
        subjectData.pMB = var.pMB;
        subjectData.pMF = var.pMF;
        %subjectData.pMTbA = var.pMTbA; % model based value transfer (based on last choice of current rated shape win/common)
        %subjectData.pMFbA = var.pMFbA; % model free value transfer (based on last choice of current rated shape win/common)
        %subjectData.ratingschange_chosens = ratingschange_chosens;

        GroupTable = [GroupTable; subjectTable];
        
        %clear subjectTable
        GroupData = [GroupData;subjectData];
        subjectRatings = subjectData.Rating;
        
        % Variables needed potentially for making graphs
        nl_ratings = table(subjectRatings, subjectData.nl_orientation);
        allratings = [allratings; nl_ratings]; % collect ratings together with nl orientation for zscoring
        %allbarv(i,:)=var.barv; % barv is the winc/common split of choices
        %allbaraa(i,:)= var.baraa; %win/common split of shape ratings, chosen not re-seen
        %allpMB(i,:)=var.pMB; % accumulare pMB from all subjects
        %allpMF(i,:)= var.pMF; % accumulate   pMF from all subjects
    
        vars = struct2table(var, 'AsArray', 1);
        grp = {getGroupLabel(group)};
        
        % performance measures
        pcorrect = nanmean(subjectData.correct);
        pstay = nanmean(subjectTable.stick); % probability of staying

        s = table(i, grp, pcorrect, pstay); % gathering data for a short table
        ss = [s, vars, dd];
        shortTable = [shortTable; ss];
        
        
    
    end % Exits the subject-specific loop
    
    GroupTable.normratings(GroupTable.nl == 1) = nanzscore(allratings.subjectRatings(allratings.Var2 == 1));
    GroupTable.normratings(GroupTable.nl == 0) = nanzscore(allratings.subjectRatings(allratings.Var2 == 0));
   
    % Gather all group data into stayprobTable
    stayprobTable = [stayprobTable; GroupTable];
    clear GroupTable
    
    %[output, ratingsTable] = AnalyseBehaviour_nl(GroupData, group, 0); % generates Daw plots for nl = 1 | nls ombined | nl = 0 ratings
    [output, ratingsTable] = AnalyseBehaviour_K1(GroupData, group, 0); % Generate plots (data, group, plot graphs(1/0))
    
    % Save each group's data into a separate file
    if options.savetables == 1
        saveFileName = sprintf("groupData_group%d.mat", group); 
        save(saveFileName, "GroupData");
    end
    
    allratings_g{group} = allratings;    % used for the smoothed plot

    data = GroupData;

    if options.optimize == 1
        fit(group) = fitWrapper(GroupData, 'LLmodelRating_K', 1, 1); %data, fucntion, EM, parallel
    else
        load('fit_new.mat')
        fit = fit(group); 
    end
    
    % Generate surrogate data
    if options.GenerateSurrData == 1
        SurrogateData = generateSurrData_K(data,fit.results.paramfit,'LLmodelRating_K'); 
        FileName = sprintf("SurrData_Group%d", group);
        save(FileName, "SurrogateData");
    else 
        load(sprintf('SurrData_Group%d.mat', group));   
    end
    
    for subject = 1:length(GroupData)
        y = [fit.results.paramfit(subject, :)];
        [LL, RPEs] = LLmodelRating_extraction1(y, data, subject); % generate RPEs based on group data

        % Make a table from the RPE data to append to stayprobTable
        subjectRPEt = struct2table(RPEs); % expand the structure into a table
        params = [params; y];
        groupRPEt = [groupRPEt; subjectRPEt]; % aggregate the RPE data within the group
        clear subjectRPEt % to prevent overwriting
    end
    
    allRPEt = [allRPEt; groupRPEt]; % aggregate all the RPEs into all one table 
    clear groupRPEt
    
    if options.fitRegression == 1   % fits a linear regression of RPEs vs ratingschange                                                                                                                                                                             
        regressions{group} = fitRegression_April(SurrogateData, data, group);
    end
end % Exits the group loops

stayprobTable = [stayprobTable, allRPEt]; % append the RPE table to stayprobTable

%% Section from LMM_table script
% stayprobTable.selected = nanzscore(stayprobTable.selected);
% stayprobTable.con_lastseen = nanzscore(stayprobTable.con_lastseen);
% stayprobTable.win_lastseen = nanzscore(stayprobTable.win_lastseen);
stayprobTable.tp = categorical(stayprobTable.tp);

TWcontrols = stayprobTable((stayprobTable.disease == 0),:);
TWpdon = stayprobTable((stayprobTable.medication == 1),:);
TWpdoff = stayprobTable((stayprobTable.medication == 0 & stayprobTable.disease == 1),:);
TWpd = stayprobTable((stayprobTable.disease == 1),:);

stats = struct();

% lmes that I want to test
model_specs = {
    % pure reward increases stay behaviour 
    'PDoffBehavLEWin', 'TWpdoff', 'stick ~ win*le + (1|subjectID)';
    'PDonBehavLEWin', 'TWpdon', 'stick ~ win*le + (1|subjectID)';
    'ControlBehavWin', 'TWcontrols', 'stick ~ win + (1|subjectID)';
    
    % model-based lmes
    'PDBehav', 'TWpd', 'stick ~ win*con*medication + (1|subjectID)';
    'ControlBehav', 'TWcontrols', 'stick ~ win*con + (1|subjectID)';
    'PDonBehav', 'TWpdon', 'stick ~ win*con + (1|subjectID)';
    'PDoffBehav', 'TWpdoff', 'stick ~ win*con + (1|subjectID)';
    'AllBehav', 'stayprobTable', 'stick ~ win*con*group + (1|subjectID)';

    % effects of clinical variables on stay behaviour/choice behaviour
    'PDBehav_pdseverity', 'TWpd', 'stick ~ win*con*medication*pd_severity + (1|subjectID)';
    'Behav_ami', 'stayprobTable', 'stick ~ win*con*group*motivation + (1|subjectID)';
    'Behav_hads', 'stayprobTable', 'stick ~ win*con*group*depression + (1|subjectID)';
    'Behav_grit', 'stayprobTable', 'stick ~ win*con*group*grit + (1|subjectID)';
    'Behav_impulsivity', 'stayprobTable', 'stick ~ win*con*group*impulsivity + (1|subjectID)';

    % all group model of rating strategy
    'PO1_ALL', 'stayprobTable', 'normratings ~ win_lastseen*con_lastseen*selected*group + (1|subjectID)';
    'PO1_ALL', 'stayprobTable', 'normratings ~ win_lastseen*con_lastseen*selected*group*nl + (1|subjectID)';

    % single group models of rating strategy (simplest)
    'PO1_HC', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected + (1|subjectID)';
    'PO1_PDon', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected + (1|subjectID)';
    'PO1_PDoff', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected + (1|subjectID)';

    % single group models of rating strategy accounting for number line effects
    'PO1_HC_nl', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected*nl + (1|subjectID)';
    'PO1_PDon_nl', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*nl + (1|subjectID)';
    'PO1_PDoff_nl', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*nl + (1|subjectID)';

    % single group models of rating strategy accounting for reaction time
    'PO1_HC_rttp', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected*RTTP + (1|subjectID)';
    'PO1_PDon_rttp', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*RTTP + (1|subjectID)';
    'PO1_PDoff_rttp', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*RTTP + (1|subjectID)';

    % effect of clinical variables on the single group rating models
    'PO1_HC_ami', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected*motivation + (1|subjectID)';
    'PO1_PDon_ami', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*motivation + (1|subjectID)';
    'PO1_PDoff_ami', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*motivation + (1|subjectID)';

    'PO1_HC_hads', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected*depression + (1|subjectID)';
    'PO1_PDon_hads', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*depression + (1|subjectID)';
    'PO1_PDoff_hads', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*depression + (1|subjectID)';
   
    'PO1_HC_grit', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected*grit + (1|subjectID)';
    'PO1_PDon_grit', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*grit + (1|subjectID)';
    'PO1_PDoff_grit', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*grit + (1|subjectID)';

    'PO1_HC_impulsivity', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected*impulsivity + (1|subjectID)';
    'PO1_PDon_impulsivity', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*impulsivity + (1|subjectID)';
    'PO1_PDoff_impulsivity', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*impulsivity + (1|subjectID)';
    
    'PO1_PDon_pdseverity', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected*pd_severity + (1|subjectID)';
    'PO1_PDoff_pdseverity', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected*pd_severity + (1|subjectID)';

    'PO1_PDoff_pdseverity', 'TWpd', 'normratings ~ win_lastseen*con_lastseen*selected*pd_severity*medication + (1|subjectID)';
    
    % added effect of trial to single group rating models
    'PO1_HC_T', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (1|subjectID)';
    'PO1_PDon_T', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (1|subjectID)';
    'PO1_PDoff_T', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (1|subjectID)';

    % added random slope of the rated shape (should it also be a fixed effect?)
    'PO1_HC_T_tpRS', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (tp|subjectID)';
    'PO1_PDon_T_tpRS', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (tp|subjectID)';
    'PO1_PDoff_T_tpRS', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (tp|subjectID)';
    
    % random effect of rated shape
    'PO1_HC_T_tpRE', 'TWcontrols', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (1|tp) + (1|subjectID)';
    'PO1_PDon_T_tpRE', 'TWpdon', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (1|tp) + (1|subjectID)';
    'PO1_PDoff_T_tpRE', 'TWpdoff', 'normratings ~ win_lastseen*con_lastseen*selected + trial + (1|tp) + (1|subjectID)';
    
    % temporal discounting factor of last seen distance 
    'PO1_HC_T_td', 'TWcontrols', 'normratings ~ selected*win_lastseen*con_lastseen*lastseendis + trial + (1|subjectID)';
    'PO1_PDon_T_td', 'TWpdon', 'normratings ~ selected*win_lastseen*con_lastseen*tp*lastseendis + trial + (1|subjectID)';
    'PO1_PDoff_T_td', 'TWpdoff', 'normratings ~ selected*win_lastseen*con_lastseen*tp*lastseendis + trial + (1|subjectID)';

     % single group models of rating strategy accounting for reaction time
    'HC_simplewin', 'TWcontrols', 'normratings ~ win*con + trial + (1|subjectID)';
    'PDon_simplewin', 'TWpdon', 'normratings ~ win*con + trial + (1|subjectID)';
    'PDoff_simplewin', 'TWpdoff', 'normratings ~ win*con + trial + (1|subjectID)';
};

% Loop through each model specification
for i = 1:size(model_specs, 1)
    model_name = model_specs{i, 1}; % get the model name from the specified model
    data_var = model_specs{i, 2}; % get the data table specified (Pd on/off, controls, all...)
    formula = model_specs{i, 3}; % get the model formula
    
   try
        % Fit the model
        d=eval(data_var);
        model = fitlme( d, formula);       
        
        % Extract statistics using a self-made function
        stats = extract_stats(model, i, stats, model_name, data_var);
        stats.fm = formula(i);
   catch ME % to catch is i've made a mistake in the formula
        fprintf('Error with model: %s\n', model_name);
        fprintf('Error message: %s\n', ME.message);
    end
end

LMMT = struct2table(stats); % convert the model stats to a table
hs = size(LMMT);

% Initialize empty tables for significant effects
significantFEs = table();
significantFEs_stick = table();
significantFEs_normratings = table();

for i = 1:hs(1)
    FEs = LMMT.FEs{i,:};
    data = LMMT.Data(i);
    formula = model_specs(i,3);
    for j = 1:length(FEs)
        if FEs(j).pValue < 0.05
            f = struct2table(FEs(j,:));
            f.Name = string(f.Name);
            f.data = string(data);
            f.AIC = LMMT.AIC(i);
            f.BIC = LMMT.BIC(i);
            f.LL = LMMT.LogLikelihood(i);
            f.formula = string(formula);
            
            significantFEs = [significantFEs; f];

            % Check if the model explains 'stick' or 'normratings' or 'RPEs'
            if contains(formula, 'stick')
                significantFEs_stick = [significantFEs_stick; f];
            elseif contains(formula, 'normratings')
                significantFEs_normratings = [significantFEs_normratings; f];
            elseif contains(formula, 'RPE')
                signifianceFEs_RPEs = [significantFEs_RPEs; f]; % significant RPE effects
            end
        else
            continue; % If any pValue in FEs is less than 0.05, add the row and break
        end
    end
end

% Convert arrays to tables if not empty
if ~isempty(significantFEs_stick)
    writetable(significantFEs_stick, 'significantFEs_stick.xlsx');
end

if ~isempty(significantFEs_normratings)
    writetable(significantFEs_normratings, 'significantFEs_normratings_new.xlsx');
end

% Write the combined table to an Excel file
if options.savetables == 1
    writetable(significantFEs, 'significantFEs.xlsx');
end



%% Pipeline for establishing the best fitting RPE
models = struct();
rpeNames = {'RPE1', 'RPE2', 'RPE3', 'RPE4', 'RPE5', 'RPE6', 'RPE7', 'RPE8', 'RPE9', 'RPE10', 'RPE11'}; % Your list of RPEs
% accumulate the linear models in models structure
for i = 1:length(rpeNames)
    formula = ['ratingschange_chosens ~ ' rpeNames{i} ' + (1|subjectID)']; % Construct a formula string
    models.(rpeNames{i}) = fitlme(stayprobTable, formula);
end

aicValues = zeros(length(rpeNames), 1);
bicValues = zeros(length(rpeNames), 1);

for i = 1:length(rpeNames)
    aicValues(i) = models.(rpeNames{i}).ModelCriterion.AIC;
    bicValues(i) = models.(rpeNames{i}).ModelCriterion.BIC;
    LL(i) = models.(rpeNames{i}).LogLikelihood;
end

% Find the model with the minimum AIC and BIC
[minAIC, idxAIC] = min(aicValues);
[minBIC, idxBIC] = min(bicValues);
[minLL, idxLL] = min(LL);

bestModelAIC = rpeNames{idxAIC};
bestModelBIC = rpeNames{idxBIC};
bestModelLL = rpeNames{idxLL};

%% determine the best fitting RPE for every group indidvidually
for i = 1:length(rpeNames)
    formula = ['ratingschange_chosens ~ ' rpeNames{i} ' + (1|subjectID)']; % Construct a formula string
    %models.(rpeNames{i}) = 
    fitlme(TWcontrols, formula)
    % to study a different group, just change the sub-table used:
    % TWcontrols, TWpdon, TWpdoff
end

% determine the best fit
aicValues = zeros(length(rpeNames), 1);
bicValues = zeros(length(rpeNames), 1);
for i = 1:length(rpeNames)
    aicValues(i) = models.(rpeNames{i}).ModelCriterion.AIC;
    bicValues(i) = models.(rpeNames{i}).ModelCriterion.BIC;
    LL(i) = models.(rpeNames{i}).LogLikelihood;
end

% Find the model with the minimum AIC and BIC
[minAIC, idxAIC] = min(aicValues);
[minBIC, idxBIC] = min(bicValues);
[minLL, idxLL] = min(LL);

bestModelAIC = rpeNames{idxAIC};
bestModelBIC = rpeNames{idxBIC};
bestModelLL = rpeNames{idxLL};