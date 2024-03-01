%% First analyse effect of reward * trial on STAY behaviour
% for on/off comparison, without PD116
klara = stayprobTable.subjectID == 1;
controls = (stayprobTable.disease == 0 & stayprobTable.subjectID ~= 35);
pdon = (stayprobTable.medication == 1 & stayprobTable.subjectID ~= 15);
pdoff = (stayprobTable.medication == 0 & stayprobTable.disease == 0 & stayprobTable.subjectID ~= 15);
pd = (stayprobTable.disease == 1);
excludedPs = (stayprobTable.subjectID == 15);
SingleSubject = stayprobTable.disease == 1 & stayprobTable.subjectID == 6; 

%KlaraBehav = fitlme(stayprobTable(klara,:), 'normratings ~ win*con*trial')

PDRatings = fitlme(stayprobTable, 'chosens ~ seen_between_chosen_ratings + (1|subjectID)')
PDoffBehav = fitlme(stayprobTable(pdoff,:), 'stick ~ win*con + (1|subjectID)')
PDonBehav = fitlme(stayprobTable(pdon,:), 'stick ~ win*con + (1|subjectID)')
HCBehav = fitlme(stayprobTable(controls,:), 'stick ~ win*con + (1|subjectID)')

PDoffRatinshchange = fitlme(stayprobTable(pdoff,:), 'ratingschange_chosens ~ win*con + trial + (1|subjectID)')
PDonRatinshchange = fitlme(stayprobTable(pdon,:), 'ratingschange_chosens ~ win*con + trial + (1|subjectID)')
ControlsRatinshchange = fitlme(stayprobTable(controls,:), 'ratingschange_chosens ~ win*con + trial + (1|subjectID)')

PDoffRatings = fitlme(stayprobTable(pdoff,:), 'normratings ~ win*con + (1|subjectID)')
PDonRatings = fitlme(stayprobTable(pdon,:), 'normratings ~ win*con + (1|subjectID)')
ControlRatings = fitlme(stayprobTable(controls,:), 'normratings ~ win*con + (1|subjectID)')


% onFirst = stayprobTable.SubjectID < 10 & pdoff;
% offFirst = stayprobTable.SubjectID >= 10 & pdoff;

%Fullmodel = fitlme(stayprobTable, 'Stick ~ 0+Medication+Disease+Reward*Consistency + (1|SubjectID)')
% pdmodel = fitlme(stayprobTable, 'Stick ~ Medication*Disease*Reward*Consistency + (1|SubjectID)')
% %outlierModel = fitlme(stayprobTable(excludedPs, :), 'Stick ~ Reward*Consistency*Medication')
% 
% ControlModel = fitlme(stayprobTable(Controls,:), 'Stick ~ Reward*Consistency + (1|SubjectID)')
% PdonModel = fitlme(stayprobTable(pdon,:), 'Stick ~ Reward*Consistency + (1|SubjectID)')
% PdoffModel = fitlme(stayprobTable(pdoff,:), 'Stick ~ Reward*Consistency + (1|SubjectID)')
% 
% Controlratings = fitlme(stayprobTable(Controls,:), 'Ratings ~ Trial*Reward*Consistency + (1|SubjectID)')
PdonNormratings = fitlme(stayprobTable(pdon,:), 'normratings ~ trial*win*con + (1|subjectID)')
Pdonratings = fitlme(stayprobTable(pdon,:), 'Ratings ~ Trial*Reward*Consistency + (1|SubjectID)')


% Pdoffratings = fitlme(stayprobTable(pdoff,:), 'Ratings ~ Trial*Reward*Consistency + (1|SubjectID)')
% 
% Fullratings = fitlme(stayprobTable(pd,:), 'Ratings ~ Medication*Trial*Reward*Consistency + (1|SubjectID)')

% % Ratinshchange models (Ratings 1-9 are discarded)
% RCmodelPDon = fitlme(stayprobTable(pdon,:), 'Ratingschange ~ Reward*Consistency + (1|SubjectID)')
% RCmodelPDoff = fitlme(stayprobTable(pdoff,:), 'Ratingschange ~ Reward*Consistency + (1|SubjectID)')
% RCmodelControls = fitlme(stayprobTable(Controls,:), 'Ratingschange ~ Reward*Consistency + (1|SubjectID)')
% RCmodelFull = fitlme(stayprobTable, 'Ratingschange ~ Reward*Consistency + (1|SubjectID)')
% 
% % ratingschange models for the current rated shape rating - prev rating,
% % grouped based on the last time it was chosen in between these two ratings
% 
% DiffRatPDon = fitlme(stayprobTable(pdon,:), 'Diffratings ~ Wons*Consis + (1|SubjectID)')
% DiffRatPDoff = fitlme(stayprobTable(pdoff,:), 'Diffratings ~ Wons*Consis + (1|SubjectID)')
% DiffRatPControls = fitlme(stayprobTable(Controls,:), 'Diffratings ~ Wons*Consis + (1|SubjectID)')
% DiffRatFull = fitlme(stayprobTable, 'Diffratings ~ Wons*Consis + (1|SubjectID)')

