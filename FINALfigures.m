%% Final Figures
tf = 16; % set the font for all titles

% Figure 2: Behaviour figure ( Behav | Ratings | Location | MBi MFi )
for group = 1:3
    data = load(sprintf('groupData_group%d', group));
    data = data.GroupData;
    ouput = AnalyseBehaviour_K1(data, group,0);
end



%% Extracting Random effect of trial 
% pd off random effects of trial
Tpdoff = fitlme(TWpdoff, 'normratings ~ win + trial + (trial|subjectID)')
[~,~,OFFstats] = randomEffects(Tpdoff);
Trial_estimates{2} = OFFstats(2:2:end,:);

% pd on random effects of trial
Tpdon = fitlme(TWpdon, 'normratings ~ win + trial + (trial|subjectID)')
[~,~,ONstats] = randomEffects(Tpdon);
Trial_estimates{1} = ONstats(2:2:end,:);

% controls random effect of trial 
Tcontrols = fitlme(TWcontrols, 'normratings ~ win + trial + (trial|subjectID)')
[~,~,Cstats] = randomEffects(Tcontrols);
Trial_estimates{3} = Cstats(2:2:end,:);

column = [];
for i = 1:3 
    column = [column; Trial_estimates{i}(:,4)];
end 
t = dataset2struct(column); % follow by copy past the structure into a new column 
tt = struct2table(t);

shortTable.Estimate = tt.Estimate;


% Figure 3: RE of Trial and Correlation between AMI Subscores & RE Trail


X = shortTable.AMITotal;
X2 = shortTable.AMIBehaviour;
X3 = shortTable.AMISocial;
X4 = shortTable.AMIEmotional;
Y = shortTable.Estimate;
Z = shortTable.grp;
ID = shortTable.i;

% Initialize gramm objects
g(1,1) = gramm('x', Z, 'y', Y, 'color', Z, 'group', Z); % Bar chart of Estimates with error bars
g(1,2) = gramm('x', X, 'y', Y, 'color', Z, 'group', Z); % AMI Total plot
g(1,3) = gramm('x', X2, 'y', Y, 'color', Z, 'group', Z); % AMI Social plot
g(1,4) = gramm('x', X3, 'y', Y, 'color', Z, 'group', Z); % AMI Social plot
g(1,5) = gramm('x', X4, 'y', Y, 'color', Z, 'group', Z); % AMI Emotional plot

% Settings for Bar chart (Estimates with error bars)
g(1,1).geom_jitter('width', 0.4); % Jitter plot with dodging
g(1,1).stat_summary('type', 'sem', 'geom', {'point', 'black_errorbar'}); % Add error bars
g(1,1).set_names('x', 'Group', 'y', 'Random Effect of Trial (x10e-3)', 'color', 'Group'); % Set axis labels
g(1,1).geom_hline('yintercept', 0, 'style', 'k--');
g(1,1).set_title('Random Effects of Trial', 'FontSize', tf)

% Settings for AMI plot
g(1,2).geom_point(); 
g(1,2).geom_hline('yintercept', 0, 'style', 'k--');  % Add a dashed horizontal line at y = 0
g(1,2).stat_glm('fullrange', true, 'geom', {'line', 'area'}); % Add summary statistics or glms
g(1,2).set_names('x', 'AMI Score', 'y', '', 'color', 'Group'); % Set axis labels
g(1,2).set_title('AMI', 'FontSize', tf); % Set plot title

% Settings for Grit plot
g(1,3).geom_point(); 
g(1,3).geom_hline('yintercept', 0, 'style', 'k--');  % Add a dashed horizontal line at y = 0
g(1,3).stat_glm('fullrange', true, 'geom', {'line', 'area'}); % Add summary statistics or glms
g(1,3).set_names('x', 'AMI Behaviour Score', 'y', '', 'color', 'Group'); % Set axis labels
g(1,3).set_title('Behaviour', 'FontSize', tf); % Set plot title

% Settings for AMI Social plot
g(1,4).geom_point(); 
g(1,4).geom_hline('yintercept', 0, 'style', 'k--');  % Add a dashed horizontal line at y = 0
g(1,4).stat_glm('fullrange', true, 'geom', {'line', 'area'}); % Add summary statistics or glms
g(1,4).set_names('x', 'AMI Social Score', 'y', '', 'color', 'Group'); % Set axis labels
g(1,4).set_title('Social', 'FontSize', tf); % Set plot title

% Settings for AMI Emotional plot
g(1,5).geom_point(); 
g(1,5).geom_hline('yintercept', 0, 'style', 'k--');  % Add a dashed horizontal line at y = 0
g(1,5).stat_glm('fullrange', true, 'geom', {'line', 'area'}); % Add summary statistics or glms
g(1,5).set_names('x', 'AMI Emotional Score', 'y', '', 'color', 'Group'); % Set axis labels
g(1,5).set_title('Emotional', 'FontSize', tf); % Set plot title
for i = 1:5
    g(1,i).no_legend()
end

% Overall settings for the combined figure
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16); % Adjust font size

% Draw the plots
f = figure('PaperOrientation','landscape','Position', [100 100 1600 400]);
g.draw();
print(gcf, 'figures/random_effects_plot', '-dpdf', '-bestfit');




% Figure 4: Ratings over time | Ratings over time split by high / low AMI (Contols and PD)
figure('PaperOrientation','landscape','Position', [100 100 1600 400]);
% Plot of ratings for all three groups
subplot(1,3,1)

conditionalPlot(TWcontrols.trial, TWcontrols.normratings, [], 'color', [1.00,0.37,0.41])
hold on
conditionalPlot(TWpdon.trial, TWpdon.normratings, [], 'color', [0.00,0.66,1.00])
conditionalPlot(TWpdoff.trial, TWpdoff.normratings, [], 'color', [0.03,0.71,0.29])
line(get(gca, 'XLim'), [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
hold off
legend({'','Controls','','PD on','','PD off'}, 'Fontsize', tf-6); legend('boxoff')
title('Ratings over trials', 'FontSize', tf); ylabel('Ratings')
xlabel('Trial', 'FontSize', tf-4)
ylabel('Ratings', 'FontSize', tf-4)
ylim([-0.5 0.5])


% High and low AMI in control
%motivated people start off lower, but for every trial motivated people
%increase ratings by .007
subplot(1,3,2)
mdl = fitlme(TWcontrols, 'normratings ~ trial*motivation*nl + (1|subjectID)');
group = TWcontrols.motivation < median(TWcontrols.motivation);
conditionalPlot(TWcontrols.trial(group), TWcontrols.normratings(group), [], 'color', [1, 0.6, 0.6])
hold on
conditionalPlot(TWcontrols.trial(~group), TWcontrols.normratings(~group), [], 'color', [0.8, 0, 0])
line(get(gca, 'XLim'), [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
hold off
xlabel('Trial', 'FontSize', tf-4)
ylabel('Ratings', 'FontSize', tf-4); ylim([-0.5, 0.5]);
legend({ "", 'Low AMI', "",'High AMI'}, 'Fontsize', tf-6); legend('boxoff')
title('Controls', 'FontSize',tf)

% Effect of high / low AMI in PD

group = TWpd.motivation < median(TWpd.motivation);
med = TWpd.medication == 1;

mdl = fitlme(TWpd, 'normratings ~ trial*motivation*medication + (1|subjectID)');
% three way interaction indicates that  apathetic patients when off reduce their ratings, but the drug prevents this
subplot(1,3,3)
conditionalPlot(TWpd.trial(group & med), TWpd.normratings(group & med), [], 'color', [0, 0, 0.5]) % dark blue, pd on
hold on
conditionalPlot(TWpd.trial(~group & med), TWpd.normratings(~group & med), [], 'color', [0.6, 0.6, 1]) % light blue, pd on
conditionalPlot(TWpd.trial(group & ~med), TWpd.normratings(group & ~med), [], 'color', [0, 0.5, 0]) % dark green, pd off
conditionalPlot(TWpd.trial(~group & ~med), TWpd.normratings(~group & ~med), [], 'color', [0.6, 1, 0.6]) % light green, pd off
line(get(gca, 'XLim'), [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
hold off
xlabel('Trial', 'FontSize', tf-4)
ylabel('Ratings', 'FontSize', tf-4); ylim([-0.5, 0.5])
legend({"",'ON: Low AMI',  "",'ON: High AMI',"", 'OFF: Low AMI', "",'OFF: High AMI'}, 'Fontsize', tf-6); legend('boxoff');
title('PD', 'FontSize',tf)
print(gcf, 'figures/ratings plots', '-dpdf', '-bestfit');

% Figure: Performance 
clear g;

% Define data 
GRP = shortTable.grp;
MF = shortTable.pMF;
MB = shortTable.pMB;
AMIT = shortTable.AMITotal;
AMIB = shortTable.AMIBehaviour;
AMIS = shortTable.AMISocial;
AMIE = shortTable.AMIEmotional;
GritTotal = shortTable.GritTotal;
SUPPSP = shortTable.SUPPSP;
HADSTotal = shortTable.HADSTotal;
pstay = shortTable.pstay;
pcorrect = shortTable.pcorrect;
pwin = shortTable.pwin;
id = shortTable.i;

% MF: Create gramm objects for the first set of plots
g(1,1) = gramm('x', GRP, 'y', pwin, 'color', GRP, 'group', GRP);
g(1,2) = gramm('x', GRP, 'y', pstay, 'color', GRP, 'group', GRP);
g(1,3) = gramm('x', GRP, 'y', pcorrect, 'color', GRP, 'group', GRP);

% Apply methods to each gramm object separately for the first set of plots

for j = 1:3
    g(1,j).geom_jitter();
     g(1,j).stat_summary('type','sem', 'geom', 'black_errorbar')
end

% Set plot labels for each subplot 
g(1,1).set_names('x', ' ', 'y', 'Performance (% wins)');
g(1,2).set_names('x', ' ', 'y', 'Performance (% correct)');
g(1,3).set_names('x', ' ', 'y', 'Stay probability');

% Remove all legends and titles
g.set_layout_options('legend', false);
g.set_title('Performance');

% Create a new figure window for the combined plots
figure('Name', 'Performance plots', 'PaperOrientation', 'landscape', 'Position', [200 400 1600 400]);

% Draw the combined plots
g.draw();
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);

%% Supplementary figures
clear g;

% Define data for the first set of plots
GRP = shortTable.grp;
ID = shortTable.i;
a = shortTable.pMF;
b = shortTable.pMB;
l = shortTable.AMITotal;
p = shortTable.AMIBehaviour;
o = shortTable.AMISocial;

% MF: Create gramm objects for the first set of plots
g(1,1) = gramm('x', GRP, 'y', a, 'color', GRP, 'subset', ID);
g(1,2) = gramm('x', GRP, 'y', b, 'color', GRP, 'subset', ID);
g(1,3) = gramm('x', GRP, 'y', l, 'color', GRP, 'subset', ID);
g(1,4) = gramm('x', GRP, 'y', p, 'color', GRP, 'subset', ID);
g(1,5) = gramm('x', GRP, 'y', o, 'color', GRP, 'subset', ID);

% Apply methods to each gramm object separately for the first set of plots
for i = 1:5
    g(1,i).geom_jitter();
    g(1,i).stat_glm('geom', 'line', 'fullrange', 1, 'disp_fit', 0);  % CI and bars
    g(1,i).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
end

% Set plot labels for each subplot in the first set of plots
g(1,1).set_names('x', 'AMI Total Score', 'y', 'Model Free Index');
g(1,2).set_names('x', 'AMI Behaviour Score', 'y', 'Model Free Index');
g(1,3).set_names('x', 'AMI Social Score', 'y', 'Model Free Index');
g(1,4).set_names('x', 'AMI Emotional Score', 'y', 'Model Free Index');

g(2,1).set_names('x', 'AMI Total Score', 'y', 'Model Based Index');
g(2,2).set_names('x', 'AMI Behaviour Score', 'y', 'Model Based Index');
g(2,3).set_names('x', 'AMI Social Score', 'y', 'Model Based Index');
g(2,4).set_names('x', 'AMI Emotional Score', 'y', 'Model Based Index');

% Define data for the second set of plots
GRP = shortTable.grp;
MF = shortTable.pMF;
MB = shortTable.pMB;
GritTotal = shortTable.GritTotal;
SUPPSP = shortTable.SUPPSP;
HADSTotal = shortTable.HADSTotal;
pstay = shortTable.pstay;
pcorrect = shortTable.pcorrect;

% MF: Create gramm objects for the second set of plots
g(3,1) = gramm('x', GritTotal, 'y', MF, 'color', GRP);
g(3,2) = gramm('x', SUPPSP, 'y', MF, 'color', GRP);
g(3,3) = gramm('x', HADSTotal, 'y', MF, 'color', GRP);
g(3,4) = gramm('x', pstay, 'y', MF, 'color', GRP);
g(3,5) = gramm('x', pcorrect, 'y', MF, 'color', GRP);

% MB: Create gramm objects for the second set of plots
g(4,1) = gramm('x', GritTotal, 'y', MB, 'color', GRP);
g(4,2) = gramm('x', SUPPSP, 'y', MB, 'color', GRP);
g(4,3) = gramm('x', HADSTotal, 'y', MB, 'color', GRP);
g(4,4) = gramm('x', pstay, 'y', MB, 'color', GRP);
g(4,5) = gramm('x', pcorrect, 'y', MB, 'color', GRP);

% Apply methods to each gramm object separately for the second set of plots
for i = 3:4
    for j = 1:5
        g(i,j).geom_point();
        g(i,j).stat_glm('geom', 'line', 'fullrange', 1, 'disp_fit', 0);  % CI and bars
        g(i,j).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
    end
end

% Set plot labels for each subplot in the second set of plots
g(3,1).set_names('x', 'Grit Score', 'y', 'Model Free Index');
g(3,2).set_names('x', 'Impulsivity', 'y', 'Model Free Index');
g(3,3).set_names('x', 'HADS Score', 'y', 'Model Free Index');
g(3,4).set_names('x', 'Probability Stay', 'y', 'Model Free Index');
g(3,5).set_names('x', 'Probability Correct', 'y', 'Model Free Index');

g(4,1).set_names('x', 'Grit Score', 'y', 'Model Based Index');
g(4,2).set_names('x', 'Impulsivity', 'y', 'Model Based Index');
g(4,3).set_names('x', 'HADS Score', 'y', 'Model Based Index');
g(4,4).set_names('x', 'Probability Stay', 'y', 'Model Based Index');
g(4,5).set_names('x', 'Probability Correct', 'y', 'Model Based Index');

% Set unified point and line options for all plots
for i = 1:4
    for j = 1:5
        g(i,j).set_point_options('base_size', 8);
        g(i,j).set_line_options('base_size', 2);
    end
end

% Remove all legends and titles
g.set_layout_options('legend', false);
g.set_title('');

% Create a new figure window for the combined plots
figure('Name', 'Combined Clinical Regressions', 'PaperOrientation', 'landscape', 'Position', [200 400 1600 1200]);

% Draw the combined plots
g.draw();
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);



