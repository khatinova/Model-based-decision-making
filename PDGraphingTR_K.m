function [pMB, pMF, pwin, barv, baraa, MBMFD, RatingsTable,...
    pMTbA, pMFbA, pMBTPOSA, pMBNEG, RTTP, baraaq, ratingschange_chosens] = PDGraphingTR_K(result, sub) 

%% PDGraphingTR 
c1=[result.A(:,1)]';
win=[result.R];
con=[result.trans];
start = [result.start];
finish = [result.finish];
tp = [result.ratingIdx]; % 1 to 5, which item probed, also referred to as ratedshape in some scripts
irrels = [result.shapeIdx]; % 
conmap=[result.trans]';
correct =[result.correct]';
startTP = [result.startTransferChoice];
endTP = [result.endTransferChoice];
normratings = [result.zRating]';
SubjectID = result.ID;
subjectID = repelem(SubjectID, length(c1), 1);
trial = (1:length(c1))';
Ratings = [result.Rating]';
LearningEffect = [result.LearningEffect]; 
medication= repelem(result.Med, length(c1), 1);
disease = repelem(result.Dis, length(c1), 1);

% PDGraphingTR continued
stick=[nan c1(1:end-1)==c1(2:end)];
nanmean(stick); %prop of trials where response is repeated
stickwin = [stick == 1 & win == 1];% trials where i stuck and won as 1;
s=stick(2:end);
stick = [s NaN]';
winf=win';
winn = win';
winf(isnan(winf))=false;
conf=con';
conf(isnan(conf))=false;

select_mbmf = [
    winf(1:end-1) & conf(1:end-1)
    winf(1:end-1) & conf(1:end-1)==0
    winf(1:end-1)==0 & conf(1:end-1)
    winf(1:end-1)==0 & conf(1:end-1)==0
    ];

nantrial=isnan(winn(1:end-1))|isnan(conmap(1:end-1));
select_mbmf = (select_mbmf .* ~nantrial)>0;


barv(1)=nanmean(s(select_mbmf(1,:)));% bar(1) is rewarded and common; 
barv(2)=nanmean(s(select_mbmf(2,:)));% bar(2) is rewarded and rare
barv(3)=nanmean(s(select_mbmf(3,:)));% bar(3) is not-rewarded and common;
barv(4)=nanmean(s(select_mbmf(4,:)));% bar(4) is not-rewarded and rare;


%% Plotting Daw graphs for each participant (now as a subplot)
    % Ratingschange subplots
    % subplot(4, 5, SubjectID); % Create a subplot in a 3x6 grid
    % hold on; 
    % for i = 1:length(diffr)
    %     yValues = diffr{i};
    % 
    %     if i == 1 || i == 3
    %         scatterColor = 'blue'; 
    %     else
    %         scatterColor = 'red'; 
    %     end
    %     % Scatter plot with x-values as the cell index (i) and y-values from the cell
    %     scatter(i * ones(size(yValues)), yValues, 'filled', 'MarkerFaceColor', scatterColor);
    % 
    %     meanValue = nanmean(yValues);
    %     plot(i, meanValue, '*', 'MarkerSize', 10, 'Color', 'black');
    % end
    % hold off;
    % 
    % % bar(rats); 
    % % Customize plot labels and title if needed
    % ylabel('Ratings Change');
    % title(['Daw Graph of Ratings Change for Subject ', num2str(SubjectID)]);
    % sgtitle(['Ratingschange on Trial Types for Group ' num2str(grp)])
    % %Replace numerical labels with x-labels
    % conditions = {'R/C', 'R/R', 'NR/C', 'NR/R'};
    % xticks(1:4);
    % xticklabels(conditions); xlim([0.5, 4.5]);

    % % Behavioral (Daw) subplots
    % subplot(4, 5, SubjectID); % Create a subplot in a 3x6 grid
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

%measures of MB and MFness of BEHAVIOUR
pMB = ((barv(1)+barv(4))-(barv(3)+barv(2)));%mb measure win common + lose uncommon - win rare + lose common
pMBPOS = (barv(1))-(barv(2)); % probability of sticking with a common vs a rare winning s1 choice
pMBNEG = (barv(4))-(barv(3)); % probability of sticking with a common vs rare losing s1 choice
pMF = ((barv(1)+barv(2))-(barv(3)+barv(4))); %mf index, total stick after win - total stick after lose
mbi = pMB - pMF; % mb - mf
pwin = mean(win); % probability of winning = mean times participant won

antic1=abs(c1-3); % 1 to 2 and 2 to 1 = opposite of c1
RTTP = endTP - startTP; %time taken for rating
RT = finish - start; % time taken at s1 choice

% initialize variables for Ratings loop
chosens = nan(length(c1), 1);
unchosens = nan(length(c1), 1);
woncs = nan(length(c1), 1);
consiscs = nan(length(c1), 1);
corr = nan(length(c1), 1);
last_chosen_distance = nan(length(c1), 1);
wonus = nan(length(c1), 1);
consisucs = nan(length(c1), 1);
last_unchosen_distance = nan(length(c1), 1);
lastseendis = nan(length(c1), 1);
ratingschange_chosens = nan(length(c1), 1);
ratingschange_unchosens = nan(length(c1), 1);
seen_between_chosen_ratings = nan(length(c1), 1);
last_rated_distance_chosens = nan(length(c1), 1);
next_rated_distance_chosens = nan(length(c1), 1);
last_rated_distance_unchosens = nan(length(c1), 1);
next_rated_distance_unchosens = nan(length(c1), 1);
seen_between_unchosen_ratings = nan(length(c1), 1);
diffratings = nan(length(c1), 1);

for i = 1:length(c1) % for each trial
    if isnan(c1(i)) % if there is no c1, make chosen shapes nan
        chosens(i) = nan;
        unchosens(i) = nan;
    else
        chosens(i) = irrels( i, c1(i) );   %find the shape for the chosen colour
        unchosens(i) = irrels(i,antic1(i)); % find the shape for the unchosen colour
    end

    % find all the distances for ratingschange_chosens and diffratings
    lastcs_ratedshape = find(chosens(1:i) == tp(i), 1, 'last'); % find the last time where the rated stimulus was chosen
    lastus_ratedshape = find(unchosens(1:i) == tp(i), 1, 'last'); %find last time at which the rated stimulus was unchosen
    last_rated_ratedshape = find(tp(1:i-1)== tp(i), 1, 'last');% when was the currently rated shape rated last?
    last_rated_chosens = find(tp(1:i-1) == chosens(i), 1, 'last');
    next_rated_chosens = find(tp(i+1:end) == chosens(i), 1, 'first') + i;
    last_rated_unchosens = find(tp(1:i-1) == unchosens(i), 1, 'last');
    next_rated_unchosens = find(tp(i+1:end) == unchosens(i), 1, 'first') + i;

    if ~isempty(lastcs_ratedshape) %if the chosen shape was chosen before this trial   
        woncs(i) = win(lastcs_ratedshape); % was the chosen shape rewarded in the last trial it was selected
        consiscs(i) = conmap(lastcs_ratedshape); % last time chosen shape was selected, was it a consistent trial?
        corr(i) = correct(lastcs_ratedshape); % did they make the correct choice (different from winning ofc) 
        last_chosen_distance(i) = i - (lastcs_ratedshape); % distance to the trial the chosen shape was last chosen
        if ~isempty(last_rated_ratedshape)
             diffratings(i) = normratings(i) - normratings(last_rated_ratedshape);
        end
    else % if the chosen shape has not been previously chosen, set the variables to Nan
        woncs(i) = nan;
        consiscs(i) = nan;        
        corr(i) = nan;
        last_chosen_distance(i) = nan;
    end 

    if ~isempty(lastus_ratedshape)
        wonus(i) = win(lastus_ratedshape); %wonus tells us of the shape not chosen was a winning shape last time it was chosen
        consisucs(i) = conmap(lastus_ratedshape); %was the last time the unchsen shape was chosen a commo trial?
        last_unchosen_distance (i) = i - lastus_ratedshape; % distance to trial i for unchosen shape 
        %diffratings_lastus = normratings(lastus_ratedshape)
    else %if the unchosen shape has not previously been chosen, set variables to NaN
        wonus(i) = nan;
        consisucs(i) = nan;
        last_unchosen_distance(i) = nan;
    end 
   lastseendis (i) = min(last_chosen_distance(i), last_unchosen_distance(i));

   if ~isempty(last_rated_chosens) && ~isempty(next_rated_chosens)
       chosen_between = find(chosens(last_rated_chosens:next_rated_chosens) == chosens(i));
       seen_between_chosen_ratings(i) = length(chosen_between) + length(find(unchosens(last_rated_unchosens:next_rated_unchosens) == unchosens(i)));
       last_rated_distance_chosens(i) = i - (last_rated_chosens);
       next_rated_distance_chosens(i) = next_rated_chosens - i;
        if length(chosen_between) > 5
          % The shape was chosen more than once between consecutive ratings, skip this trial
           ratingschange_chosens(i) = NaN;
        else 
           ratingschange_chosens(i) = normratings(next_rated_chosens) - normratings(last_rated_chosens);
        end
   end
   if ~isempty(last_rated_unchosens) && ~isempty(next_rated_unchosens)
       last_rated_distance_unchosens(i) = i - (last_rated_unchosens);
       next_rated_distance_unchosens(i) = next_rated_unchosens - i;
       unchosen_between = find(unchosens(last_rated_unchosens:next_rated_unchosens) == unchosens(i));
       %seen_between_unchosen_ratings(i) = unchosen_between + find(chosens(last_rated_chosens:next_rated_chosens) == chosens(i));
       ratingschange_unchosens(i) = normratings(next_rated_unchosens) - normratings(last_rated_unchosens);
   end

end
RatingsTable = table();

RatingsTable = table(subjectID, trial, medication, disease, stick, win, con, chosens, unchosens, tp, normratings, ratingschange_chosens, ratingschange_unchosens,  ...
    diffratings, woncs, consiscs, wonus, consisucs, last_rated_distance_chosens, next_rated_distance_chosens, last_chosen_distance, last_unchosen_distance, ...
    last_rated_distance_unchosens, next_rated_distance_unchosens, seen_between_chosen_ratings);

consistentcs=(consiscs==1); % was the rated shape last chosen on a common trial?
inconsistentcs=(consiscs==0);
consistentucs=(consisucs==1); % was the rated shape unchosen on a common trial?
inconsistentucs=(consisucs==0);

selectcs = last_chosen_distance<10; %only consider ratings for shapes chosen x trials ago 
selectucs = last_unchosen_distance<10; %only give ratings for unchosen shapes unchosen x trials ago 
chosennotreseen = last_chosen_distance < last_unchosen_distance; %gives ratings for chosen shapes that havent been shown between chosen - rating 
unchosennotreseen = last_chosen_distance > last_unchosen_distance; %gives on ratings for unchosen shapes that havent been re-shown before being rated
selandnotreseen = selectcs & chosennotreseen; % ratings of chosen shapes less than X trials ago and not re-shown
unselandnotreseen = selectucs & unchosennotreseen; % ratings of unchosen shapes less than X trials ago and not re-shown

baraa(1)= nanmean(normratings(woncs==1&consistentcs&(selandnotreseen)));
baraa(2)= nanmean(normratings(woncs==1&inconsistentcs&(selandnotreseen)));
baraa(3)= nanmean(normratings(woncs==0&consistentcs&(selandnotreseen)));
baraa(4)= nanmean(normratings(woncs==0&inconsistentcs&(selandnotreseen)));


% calculate diffratings also for this 
ucbara(1)= nanmean(normratings(wonus==1&consistentucs&(unselandnotreseen)));
ucbara(2)= nanmean(normratings(wonus==1&inconsistentucs&(unselandnotreseen)));
ucbara(3)= nanmean(normratings(wonus==0&consistentucs&(unselandnotreseen)));
ucbara(4)= nanmean(normratings(wonus==0&inconsistentucs&(unselandnotreseen)));

currRating(1)= nanmean(normratings(winf==1&conf==1));
currRating(2)= nanmean(normratings(winf==1&conf==0));
currRating(3)= nanmean(normratings(winf==0&conf==1));
currRating(4)= nanmean(normratings(winf==0&conf==0));

    % Behavioral subplots
    % subplot(4, 5, sub); % Create a subplot in a 3x6 grid
    % bar(currRating); 
    % % Customize plot labels and title if needed
    % ylabel('Normrating Trial(i)');
    % title(['Subject ', num2str(sub)]);
    % sgtitle('Ratings(i) binned by reward/common(i) - Controls + KH')
    % % Replace numerical labels with x-labels
    % conditions = {'R/C', 'R/R', 'NR/C', 'NR/R'};
    % xticks(1:4);
    % xticklabels(conditions);

pMTbA = ((baraa(1)+baraa(4))-(baraa(3)+baraa(2))); %prob of MBT to bara (general conis in conis graphh)
pMFbA = ((baraa(1)+baraa(2))-(baraa(3)+baraa(4))); % prob of MFT to bara 
pMBTPOSA = ((baraa(1))-(baraa(2))); %MBT for positive (wins) only 
pMBTNEGA = ((baraa(4))-(baraa(3))); %MBT for neg only 
% figure
% bar([pMTbA, pMFbA, pMBTPOSA, pMBTNEGA] )
% ylabel('probability');
% xticklabels({'pMTbA', 'pMFbA', 'pMBT positive', 'pMBT negative'});
% xticks(1:4)
% title('pMB transfers');


%% How ratings decay with time (temporal discounting...?)
% How many trials ago was the currently rated shape last chosen
select1 = last_chosen_distance==2;
select2 = last_chosen_distance==3;
select3 = last_chosen_distance==4;
select4 = last_chosen_distance==5;
select5 = last_chosen_distance==6;
select6 = last_chosen_distance==7;
select7 = last_chosen_distance==8;
select8 = last_chosen_distance==9;

baraaq(1)= nanmean(normratings(woncs==1&consistentcs&(select1)&chosennotreseen));
baraaq(2)= nanmean(normratings(woncs==1&inconsistentcs&(select1)&chosennotreseen));
baraaq(3)= nanmean(normratings(woncs==0&consistentcs&(select1)&chosennotreseen));
baraaq(4)= nanmean(normratings(woncs==0&inconsistentcs&(select1)&chosennotreseen));
baraaq(5)= nanmean(normratings(woncs==1&consistentcs&(select2)&chosennotreseen));
baraaq(6)= nanmean(normratings(woncs==1&inconsistentcs&(select2)&chosennotreseen));
baraaq(7)= nanmean(normratings(woncs==0&consistentcs&(select2)&chosennotreseen));
baraaq(8)= nanmean(normratings(woncs==0&inconsistentcs&(select2)&chosennotreseen));
baraaq(9)= nanmean(normratings(woncs==1&consistentcs&(select3)&chosennotreseen));
baraaq(10)= nanmean(normratings(woncs==1&inconsistentcs&(select3)&chosennotreseen));
baraaq(11)= nanmean(normratings(woncs==0&consistentcs&(select3)&chosennotreseen));
baraaq(12)= nanmean(normratings(woncs==0&inconsistentcs&(select3)&chosennotreseen));
baraaq(13)= nanmean(normratings(woncs==1&consistentcs&(select4)&chosennotreseen));
baraaq(14)= nanmean(normratings(woncs==1&inconsistentcs&(select4)&chosennotreseen));
baraaq(15)= nanmean(normratings(woncs==0&consistentcs&(select4)&chosennotreseen));
baraaq(16)= nanmean(normratings(woncs==0&inconsistentcs&(select4)&chosennotreseen));
baraaq(17)= nanmean(normratings(woncs==1&consistentcs&(select5)&chosennotreseen));
baraaq(18)= nanmean(normratings(woncs==1&inconsistentcs&(select5)&chosennotreseen));
baraaq(19)= nanmean(normratings(woncs==0&consistentcs&(select5)&chosennotreseen));
baraaq(20)= nanmean(normratings(woncs==0&inconsistentcs&(select5)&chosennotreseen));
baraaq(21)= nanmean(normratings(woncs==1&consistentcs&(select6)&chosennotreseen));
baraaq(22)= nanmean(normratings(woncs==1&inconsistentcs&(select6)&chosennotreseen));
baraaq(23)= nanmean(normratings(woncs==0&consistentcs&(select6)&chosennotreseen));
baraaq(24)= nanmean(normratings(woncs==0&inconsistentcs&(select6)&chosennotreseen));
baraaq(25)= nanmean(normratings(woncs==1&consistentcs&(select7)&chosennotreseen));
baraaq(26)= nanmean(normratings(woncs==1&inconsistentcs&(select7)&chosennotreseen));
baraaq(27)= nanmean(normratings(woncs==0&consistentcs&(select7)&chosennotreseen));
baraaq(28)= nanmean(normratings(woncs==0&inconsistentcs&(select7)&chosennotreseen));
baraaq(29)= nanmean(normratings(woncs==1&consistentcs&(select8)&chosennotreseen));
baraaq(30)= nanmean(normratings(woncs==1&inconsistentcs&(select8)&chosennotreseen));
baraaq(31)= nanmean(normratings(woncs==0&consistentcs&(select8)&chosennotreseen));
baraaq(32)= nanmean(normratings(woncs==0&inconsistentcs&(select8)&chosennotreseen));
errbaraaq = baraaq - mean(baraaq);

% MB vs MF decay (MB first MF last)
MBMFD(1)=((baraaq(1)+baraaq(4))-(baraaq(3)+baraaq(2)));
MBMFD(2)=((baraaq(5)+baraaq(8))-(baraaq(7)+baraaq(6)));
MBMFD(3)=((baraaq(9)+baraaq(12))-(baraaq(11)+baraaq(10)));
MBMFD(4)=((baraaq(13)+baraaq(16))-(baraaq(15)+baraaq(14)));
MBMFD(5)=((baraaq(17)+baraaq(20))-(baraaq(19)+baraaq(18)));
MBMFD(6)=((baraaq(21)+baraaq(24))-(baraaq(23)+baraaq(22)));
MBMFD(7)=((baraaq(25)+baraaq(28))-(baraaq(27)+baraaq(26)));
MBMFD(8)=((baraaq(29)+baraaq(32))-(baraaq(31)+baraaq(30)));
MBMFD(10)=((baraaq(1)+baraaq(2))-(baraaq(3)+baraaq(4)));
MBMFD(11)=((baraaq(5)+baraaq(6))-(baraaq(7)+baraaq(8)));
MBMFD(12)=((baraaq(8)+baraaq(10))-(baraaq(11)+baraaq(12)));
MBMFD(13)=((baraaq(13)+baraaq(14))-(baraaq(15)+baraaq(16)));
MBMFD(14)=((baraaq(17)+baraaq(18))-(baraaq(19)+baraaq(20)));
MBMFD(15)=((baraaq(21)+baraaq(22))-(baraaq(23)+baraaq(24)));
MBMFD(16)=((baraaq(25)+baraaq(26))-(baraaq(27)+baraaq(28)));
MBMFD(17)=((baraaq(29)+baraaq(30))-(baraaq(31)+baraaq(32)));

% Adding cutoff of 5 as "proximate" chosen trials, anything over as "distal"
select1 = last_chosen_distance<5;
select2 = last_chosen_distance>4;

baraaqa(1)= nanmean(normratings(woncs==1&consistentcs&(select1)&chosennotreseen));
baraaqa(2)= nanmean(normratings(woncs==1&inconsistentcs&(select1)&chosennotreseen));
baraaqa(3)= nanmean(normratings(woncs==0&consistentcs&(select1)&chosennotreseen));
baraaqa(4)= nanmean(normratings(woncs==0&inconsistentcs&(select1)&chosennotreseen));
baraaqa(5)= nanmean(normratings(woncs==1&consistentcs&(select2)&chosennotreseen));
baraaqa(6)= nanmean(normratings(woncs==1&inconsistentcs&(select2)&chosennotreseen));
baraaqa(7)= nanmean(normratings(woncs==0&consistentcs&(select2)&chosennotreseen));
baraaqa(8)= nanmean(normratings(woncs==0&inconsistentcs&(select2)&chosennotreseen));
errbaraaqa = baraaqa - mean(baraaqa);

end
