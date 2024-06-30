function  [var, RatingsTable, ratingschange_chosens] = PDGraphingTR_K1(result, group) 
 % Adapted to work for the MD young healthy controls when loaded in as pure data files

%% PDGraphingTR 
c1=[result.A(:,1)];
win=[result.R]; % reward (1) or non-reward (0)
con=[result.trans]; % same as conmap, but transposed    
start = [result.start]; % start s1
finish = [result.finish]; % end s1
tp = [result.ratingIdx]; % 1 to 5, which item probed, also referred to as ratedshape in some scripts
irrels = [result.shapeIdx]; % which s1 shapes were presented (red | blue)
conmap=[result.trans]; % common (1) or rare (0)
correct =[result.correct]; % should the s2 choice have led to a reward?
startTP = [result.startTransferChoice]; % start of rating
endTP = [result.endTransferChoice]; % end of rating
normratings = [result.zRating]; % normalised ratings
SubjectID = result.ID;
subjectID = repelem(SubjectID, length(c1), 1); % subjectID for subjectTable
trial = (1:length(c1))';
ratings = [result.Rating]; % normratings

RT = finish - start; idx = true(length(c1),1); % calculate s1 choice time and remove times over 4s
if RT>4000
    idx = false;
end
RTTP = (endTP - startTP); % time taken for rating

% Variables for the table
medication = repelem(result.Med, length(c1), 1);
disease = repelem(result.Dis, length(c1), 1);
age = repelem(result.Age, length(c1), 1);
pd_severity = repelem(result.UPDRST, length(c1), 1);
motivation = repelem(result.AMI, length(c1), 1);
depression = repelem(result.HADS, length(c1), 1);
grit = repelem(result.Grit, length(c1), 1);
impulsivity = repelem(result.Impulsivity, length(c1), 1);
le = repelem(result.LearningEffect, length(c1), 1);
nl = result.nl_orientation;
Screen = result.ScreenSize;
grp = repelem(cellstr(getGroupLabel(group)),192,1);
version = categorical(result.version);
rating_coord = result.adj_coord;

% PDGraphingTR continued
stick=[nan; c1(1:end-1)==c1(2:end)]; % did they repeat ths s1 choice?
stickwin = stick == 1 & win == 1;% trials where i stuck and won as 1;
s=stick(2:end);
stick(idx == 1) = [s; NaN]; % set stick behaviour 2:end as (1:end-1 + NaN)
var.stick = stick;
win(win == 0) = -1;
win(isnan(win))=false;
con(con == 0) = -1;
con(isnan(con))=false;

select_mbmf = [
    (win(1:end-1) == 1) & (con(1:end-1) == 1), ...
    (win(1:end-1) == 1) & (con(1:end-1) == -1), ...
    (win(1:end-1) == -1) & (con(1:end-1) == 1), ...
    (win(1:end-1) == -1) & (con(1:end-1) == -1)
]; % clasifies trials 

% stay probability
barv(1)=nanmean(stick(select_mbmf(:,1)));% bar(1) is rewarded and common; 
barv(2)=nanmean(stick(select_mbmf(:,2)));% bar(2) is rewarded and rare
barv(3)=nanmean(stick(select_mbmf(:,3)));% bar(3) is not-rewarded and common;
barv(4)=nanmean(stick(select_mbmf(:,4)));% bar(4) is not-rewarded and rare;
var.barv = barv;

%% Plotting Daw graphs for each participant (now as a subplot)
    % Behavioral (Daw) subplots
    % subplot(3, 10, SubjectID); % Create a subplot in a 3x6 grid
    % barColors = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
    % bar(barv); 
    % % Customize plot labels and title if needed
    % ylabel('Stay Behaviour');
    % title(['Subject ', num2str(SubjectID)]);
    % sgtitle('Stay behaviour Daw graphs for Controls ')
    % % Replace numerical labels with x-labels
    % colormap(barColors);
    % conditions = {'R/C', 'R/R', 'NR/C', 'NR/R'};
    % xticks(1:4);
    % xticklabels(conditions);
    % ylim([0.4 1.0]);

%measures of MB and MFness of BEHAVIOUR
var.pMB = ((barv(1)+barv(4))-(barv(3)+barv(2)));%mb measure win common + lose uncommon - win rare + lose common
var.pMBPOS = (barv(1))-(barv(2)); % probability of sticking with a common vs a rare winning s1 choice
var.pMBNEG = (barv(4))-(barv(3)); % probability of sticking with a common vs rare losing s1 choice
var.pMF = ((barv(1)+barv(2))-(barv(3)+barv(4))); %mf index, total stick after win - total stick after lose
var.mbi = var.pMB - var.pMF; % mb - mf
var.pwin = (sum(win == 1))/length(c1); % probability of winning = mean times participant won

antic1=abs(c1-3); % 1 to 2 and 2 to 1 = opposite of c1


% initialize variables for Ratings loop
chosens = nan(length(c1), 1);
unchosens = nan(length(c1), 1);
ratingschange_chosens = nan(length(c1), 1);
lastseentrial2 = nan(length(c1), 1);
lastseentrial = nan(length(c1), 1);
win_lastseen = nan(length(c1), 1);
con_lastseen = nan(length(c1), 1);
selected = nan(length(c1), 1);
win_lastseen2 = nan(length(c1), 1);
con_lastseen2 = nan(length(c1), 1);
selected2 = nan(length(c1), 1);
trialdiff = nan(length(c1), 1);
rtc = nan(length(c1), 1);
nrtRating = nan(length(c1), 1);

%% Ratingschange function
for i = 1:length(c1) % for each trial
    if c1(i) == 0 || isnan(c1(i)) % if there is no c1, make (un)chosen shapes nan
        chosens(i) = nan;
        unchosens(i) = nan;
    else
        chosens(i) = irrels( i, c1(i) );   %find the shape for the chosen colour
        unchosens(i) = irrels(i,antic1(i)); % find the shape for the unchosen colour
    end

    lst = find(chosens(1:i) == tp(i) | unchosens(1:i) == tp(i), 1, 'last');
    
    % needs to be fine tuned
    if ~isempty(lst)
        lastseentrial(i) = lst;
        
        if chosens(lastseentrial(i)) == tp(i)
            selected(i) = 1;
            
        elseif unchosens(lastseentrial(i)) == tp(i)
            selected(i) = -1;
        else
            selected(i) = 0; % This case should not happen if the find condition was true
        end
        
        win_lastseen(i) = win(lastseentrial(i)); 
        % last time it was seen, did it lead to a win? (regardless of it if was chosen or not)
        con_lastseen(i) = con(lastseentrial(i));
        % last time it was seen, did it lead to a consistent trial? (regardless of it if was chosen or not)
    end
    
    lst2 = find(chosens(1:i) == tp(i) | unchosens(1:i) == tp(i), 2, 'last');
    % looks at the trial seen before the last trial seen
    if ~isempty(lst2)
        lastseentrial2(i) = lst2(1);

        if chosens(lastseentrial2(i)) == tp(i)
            selected2(i) = 1;
            
        elseif unchosens(lastseentrial2(i)) == tp(i)
            selected2(i) = -1;
        else
            selected2(i) = 0; % This case should not happen if the find condition was true
        end
        
        win_lastseen2(i) = win(lastseentrial2(i)); 
        % last time it was seen, did it lead to a win? (regardless of it if was chosen or not)
        con_lastseen2(i) = con(lastseentrial2(i));
        % last time it was seen, did it lead to a consistent trial? (regardless of it if was chosen or not)
    end

    % For the current rated shape, what was the difference between the last rating assigned and the current rating
    lrt = find(tp(1:i-1) == tp(i), 1, 'last');
    if ~isempty(lrt)
        trialdiff(i) = normratings(i) - normratings(lrt);
    end
end

for i = 1:length(c1)
        % does halo effect persist for the next time the rated shape on a win trial is rated?
    nrt_halo = find(tp(i) == tp(i+1:end), 1, 'first');

    if ~isempty(nrt_halo) 
        if isempty(find(chosens(i) == chosens(nrt_halo), 1))
            nrtRating(i) = normratings(nrt_halo);% What was rating on the next 'halo' effect rated shape?
            rtc(i) = nrtRating(i)  - normratings(i); % what was the difference between the 1st and next?
        end
    end
end

RatingsTable = table();
RatingsTable = table(subjectID, grp, Screen, trial, le, version, age, medication, disease, nl, pd_severity, motivation, depression, grit, ...
    impulsivity, stick, win, con, chosens,  unchosens, tp, trialdiff, normratings, nrtRating, rtc, rating_coord, RT, RTTP, lastseentrial, selected, win_lastseen, con_lastseen,...
    lastseentrial2, selected2, win_lastseen2, con_lastseen2);

ws(1) = nanmean(normratings((win_lastseen==1) & (selected==1)));
ws(2) = nanmean(normratings((win_lastseen==-1) & (selected==1)));
ws(3) = nanmean(normratings((win_lastseen==1) & (selected==-1)));
ws(4) = nanmean(normratings((win_lastseen==-1) & (selected==-1)));

% Calculate standard errors for the means
se(1) = nanstd(normratings((win_lastseen==1) & (selected==1))) / sqrt(sum((win_lastseen==1) & (selected==1)));
se(2) = nanstd(normratings((win_lastseen==-1) & (selected==1))) / sqrt(sum((win_lastseen==-1) & (selected==1)));
se(3) = nanstd(normratings((win_lastseen==1) & (selected==-1))) / sqrt(sum((win_lastseen==1) & (selected==-1)));
se(4) = nanstd(normratings((win_lastseen==-1) & (selected==-1))) / sqrt(sum((win_lastseen==-1) & (selected==-1)));

sgtitle('Ratings based on last seen & NOT selected', 'Fontsize', 24)
subplot(4, 5, SubjectID); % Create a subplot in a 3x6 grid
barColors = [1 0 0; 0 0 1];
% Create the bar plot
b = bar(ws(3:4)); 
hold on;
x = 1:length(ws(3:4)); % X locations of the bars
errorbar(x, ws(3:4), se(3:4), 'k', 'linestyle', 'none', 'LineWidth', 2); 

if SubjectID == 1 || SubjectID ==6 || SubjectID ==11 || SubjectID == 16
    ylabel('Zscored ratings');
end

% Replace numerical labels 
title(['sub ' num2str(SubjectID)], 'Fontsize', 12)
colormap(barColors);
conditions = {'Win', 'Non-Win'};
xticks(1:2);
xticklabels(conditions);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 24);
ylim([-0.7 0.7])

% subplot(1, 2, 2); % Create a subplot in a 3x6 grid
% barColors = [0 0 1; 1 0 0];
% bar(ws(3:4)); 
% % Customize plot labels and title if needed
% ylabel('Zscored ratings');
% title('Last Seen Not Selected');
% % Replace numerical labels with x-labels
% colormap(barColors);
% conditions = {'Win', 'Non Win'};
% xticks(1:2);
% xticklabels(conditions);

end
