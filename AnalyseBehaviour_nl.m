function [output, ratingsTable] = AnalyseBehaviour_nl(data, group,toplot)
% Inputs:
 %toplot = 1; when data should be plotted, else 0

% groups: [1 = pd off,2 = pd on,3 = Old_HC_controls]
% Output:
% output   = structure
%     .StayData     = stay probabilities
%     .StaySim      = stay probabilities simulated data
%     .Rating       = ratings
%     .Attention    = attention effect on stay behaviour
%     .Location     = location effect on stay behaviour
%     .MB           = model-based control
%     .MF           = model-free control
%     .MBsim        = model-based control simulated data
%     .MFsim        = model-free control simulated data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = getGroupLabel(group);
N = data.Nch;


nsub        = length(data);
nRepeat     = length(data(1).A);
rr          = nan(nsub,4);
rr1          = nan(nsub,4);
rr0          = nan(nsub,4);
rrs         = nan(nsub,4,nRepeat);
missed      = nan(nsub,1);
AvRating    = nan(nsub,4);
ll          = nan(nsub,4);
ratingsTable = table();

for sub = 1:nsub

    %% Analyse Real Data
    a       = data(sub).A;
    missed(sub,1) = sum(isnan(sum(a,2)));
    idx     = 1:192;
    a       = a(idx,:);
    s       = data(sub).S(idx,:);
    r       = data(sub).R(idx);
    nl      = data(sub).nl_orientation(idx);
    trans   = data(sub).trans(idx);
    stay1    = a(2:end,1)==a(1:end-1,1);
    rr(sub,1) = nanmean(stay1 (trans(1:end-1)==1 & r(1:end-1)==1)); % CR
    rr(sub,2) = nanmean(stay1 (trans(1:end-1)==0 & r(1:end-1)==1)); % RR
    rr(sub,3) = nanmean(stay1 (trans(1:end-1)==1 & r(1:end-1)==0)); % CU
    rr(sub,4) = nanmean(stay1 (trans(1:end-1)==0 & r(1:end-1)==0)); % RU
    


    %% Analyse Rating and Attention effect
    irrelshape      = data(sub).shapeIdx(idx,:);
    ratingIdx       = data(sub).ratingIdx(idx);
    c1              = data(sub).A(idx,1);
    win             = data(sub).R(idx);
    stay            = c1(2:end,1)==c1(1:end-1,1);
    correct         = data(sub).correct(idx);
    u1              = 3 - c1;
    ratings         = data(sub).Rating(idx);

    chosenShapes    = nan(N,1);
    unchosenShapes  = nan(N,1);

    for t = 1:N
        if ~isnan(c1(t))
            chosenShapes(t) = irrelshape(t,c1(t));
            unchosenShapes(t) = irrelshape(t,u1(t));
        end
    end

    woncs   = nan(N,1);    % last time won chosen shape
    wonus   = nan(N,1);    % last time won unchosen shape 
    commonC = nan(N,1);    % transition chosen shape
    commonU = nan(N,1);    % transition unchosen shape
    last_chosen_distance    = nan(N,1);    % distance between probes chosen
    last_unchosen_distance  = nan(N,1);    % distance between probes unchosen
    correctChoice           = nan(N,1);    % choice that leads to best outcome

    for i = 1:length(c1)
        % find the last time chosen stimulus was probed
        lastcs = find(chosenShapes(1:i) == ratingIdx(i), 1, 'last'); 
        % find last time the unchosen shape was choice 1 
        lastus = find(unchosenShapes(1:i) == ratingIdx(i), 1, 'last');

        % for chosen shape
        if ~isempty(lastcs)
            woncs(i)        = win(lastcs);
            commonC(i)      = trans(lastcs);
            correctChoice(i) = correct(lastcs);
            last_chosen_distance(i) = i - lastcs;
        end

        % for unchosen shape
        if ~isempty(lastus)
            wonus(i)        = win(lastus);
            commonU(i)      = trans(lastus);
            last_unchosen_distance (i) = i - lastus;
        end

    end

    normratings = data(sub).zRating;
    
    % Compute average rating for each transition and reward
    AvRating(sub,1)= nanmean(normratings(commonC == 1 & woncs == 1 ));
    AvRating(sub,2)= nanmean(normratings(commonC == 0 & woncs == 1 ));
    AvRating(sub,3)= nanmean(normratings(commonC == 1 & woncs == 0 ));
    AvRating(sub,4)= nanmean(normratings(commonC == 0 & woncs == 0 ));

    rr1(sub,1) = nanmean(normratings(commonC == 1 & woncs == 1 & nl(1) == 1)); % CR
    rr1(sub,2) = nanmean(normratings(commonC == 0 & woncs == 1 & nl(1) == 1)); % RR
    rr1(sub,3) = nanmean(normratings(commonC == 1 & woncs == 0 & nl(1) == 1)); % CU
    rr1(sub,4) = nanmean(normratings(commonC == 0 & woncs == 0 & nl(1) == 1)); % RU

    rr0(sub,1) = nanmean(normratings(commonC == 1 & woncs == 1 & nl(1) == 0)); % CR
    rr0(sub,2) = nanmean(normratings(commonC == 0 & woncs == 1 & nl(1) == 0)); % RR
    rr0(sub,3) = nanmean(normratings(commonC == 1 & woncs == 0 & nl(1) == 0)); % CU
    rr0(sub,4) = nanmean(normratings(commonC == 0 & woncs == 0 & nl(1) == 0)); % RU   
 
    %% Analyse location effects
    location        = data(sub).LRchoice(idx);
    locationStay    = location(2:end,1) == location(1:end-1,1);
    
    ll(sub,1) = mean(locationStay (trans(1:end-1)==1 & r(1:end-1)==1)); % CR
    ll(sub,2) = mean(locationStay (trans(1:end-1)==0 & r(1:end-1)==1)); % RR
    ll(sub,3) = mean(locationStay (trans(1:end-1)==1 & r(1:end-1)==0)); % CU
    ll(sub,4) = mean(locationStay (trans(1:end-1)==0 & r(1:end-1)==0)); % RU

    id = repelem(sub, 192,1);
    srt = table(id, normratings, woncs, commonC, wonus, commonU, last_chosen_distance, last_unchosen_distance);
    ratingsTable = [ratingsTable; srt];
end

% Compute the model-based and model-free values
mb = (rr(:,1) - rr(:,2)) - (rr(:,3) - rr(:,4));
mf = (rr(:,1) + rr(:,2)) - (rr(:,3) + rr(:,4));
mb_rating = (AvRating(:,1) - AvRating(:,2)) - (AvRating(:,3) - AvRating(:,4));
mf_rating = (AvRating(:,1) + AvRating(:,2)) - (AvRating(:,3) + AvRating(:,4));


stayData = rr1;
stayLocation = rr0;

%% Plot the data
if toplot   % plot the data
%% Stay behaviour real data
    f = figure('Name', 'Daw Plots', 'PaperOrientation','landscape', 'Position', [200 400 600 400]);
    sgtitle([getGroupLabel(group) ' nl = 1 | Ratings | nl = 0'])
    axes1 = subplot(1,3,1);
    hold(axes1,'on');

    y = [nanmean(stayData(:,1:2)); nanmean(stayData(:,3:4))];
    stdev = [nanstd(stayData(:,1:2)); nanstd(stayData(:,3:4))]/sqrt(nsub);

    ngroups = size(y, 1);
    nbars = size(y, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));

    ha1 = bar(y, 'Parent',axes1,'BarWidth',1, 'LineWidth', 2);
    set(ha1(1),'DisplayName','Common','FaceColor','b', 'FaceAlpha', 1);
    set(ha1(2),'DisplayName','Rare','FaceColor','r', 'FaceAlpha', 1);
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, y(:,i), stdev(:,i),'k', 'LineStyle', 'none', 'LineWidth', 1);
    end
    box(axes1,'on');
    set(axes1,'XTick',[1 2],'XTickLabel',{'Rewarded','Unrewarded'}, 'LineWidth', 1);
    set(axes1,'Ylim',[-0.2 0.3]);
    %legend({'Common' 'Rare'})
    ylabel('rating');
    title([getGroupLabel(group)]);
 
%% Rating data
    axes2 = subplot(1,3,2);
    hold(axes2,'on');

    y = [mean(AvRating(:,1:2)); mean(AvRating(:,3:4))];
    stdev = [std(AvRating(:,1:2)); std(AvRating(:,3:4))]/sqrt(nsub);

    ngroups = size(y, 1);
    nbars = size(y, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));

    ha1 = bar(y, 'Parent',axes2,'BarWidth',1, 'LineWidth', 2);
    set(ha1(1),'DisplayName','Common','FaceColor','b', 'FaceAlpha', 1);
    set(ha1(2),'DisplayName','Rare','FaceColor','r', 'FaceAlpha', 1);
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, y(:,i), stdev(:,i),'k', 'LineStyle', 'none', 'LineWidth', 1);
    end
    box(axes2,'on');
    set(axes2,'XTick',[1 2],'XTickLabel',{'Rewarded','Unrewarded'}, 'LineWidth', 1);
    %legend({'Common' 'Rare'})
    set(axes2,'Ylim',[-0.15 0.15]);
    ylabel('Rating');
    title([getGroupLabel(group)]);

    
%% Stay behaviour Location
    axes3 = subplot(1,3,3);
    hold(axes3,'on');

    y = [nanmean(stayLocation(:,1:2)); nanmean(stayLocation(:,3:4))];
    stdev = [nanstd(stayLocation(:,1:2)); nanstd(stayLocation(:,3:4))]/sqrt(nsub);

    ngroups = size(y, 1);
    nbars = size(y, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));

    ha1 = bar(y, 'Parent',axes3,'BarWidth',1, 'LineWidth', 3);
    set(ha1(1),'DisplayName','Common','FaceColor','b', 'FaceAlpha', 1);
    set(ha1(2),'DisplayName','Rare','FaceColor','r', 'FaceAlpha', 1);
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, y(:,i), stdev(:,i),'k', 'LineStyle', 'none', 'LineWidth', 1);
    end
    box(axes3,'on');
    set(axes3,'XTick',[1 2],'XTickLabel',{'Rewarded','Unrewarded'}, 'LineWidth', 1);
    set(axes3,'Ylim',[-0.2 0.3]);
    ylabel('Stay Probability');
    title('nl = 0');
    %legend({'Common' 'Rare'})

    set(f,'units','centimeter','position',[1,5,35,7])

    exportgraphics(gcf,['ratings daw by nl ', getGroupLabel(group), '.pdf'])

end

%% Store relevant data
output.stayData     = stayData;
output.Rating       = AvRating;
output.Location     = stayLocation;
output.MB           = mb;
output.MF           = mf;
output.MBrating     = mb_rating;
output.MFrating     = mf_rating;
output.missed       = missed;

end
