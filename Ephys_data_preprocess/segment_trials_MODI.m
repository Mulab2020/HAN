function [trialData,idx] = segment_trials_MODI(data, bouts)
% data.CurrObsDur as segment criteria
disp('Segmenting...');
trial_num = data.trial_num;

trial_num(trial_num==0) = [];
idx = zeros(1,trial_num(end)-trial_num(1)+1);
trialData = cell(1,trial_num(end)-trial_num(1)+1);
for i = trial_num(1):data.trial_num(end)

    trialIdx = find(data.trial_num == i);
    trialDatum.trial_start = trialIdx(1);
    trialDatum.trial_end = trialIdx(end);
    trialDatum.trialType = data.stimParam2(floor((trialDatum.trial_start+trialDatum.trial_end)/2)); % data.CurrObsDur or stimParam2 depends
    idx(i-trial_num(1)+1) = trialDatum.trialType;
    trialDatum.bout_start = bouts.bout_start(bouts.bout_start>=trialDatum.trial_start & bouts.bout_start<= trialDatum.trial_end);
    trialDatum.bout_end = bouts.bout_end(bouts.bout_start>=trialDatum.trial_start & bouts.bout_start<= trialDatum.trial_end);
    trialDatum.bout_power = bouts.bout_power(:,bouts.bout_start>=trialDatum.trial_start & bouts.bout_start<= trialDatum.trial_end);
    trialDatum.xpos = data.fishx(trialDatum.trial_start:trialDatum.trial_end);
    trialDatum.xpos = trialDatum.xpos-trialDatum.xpos(1);
    trialDatum.ypos = data.fishy(trialDatum.trial_start:trialDatum.trial_end);
    trialDatum.ypos = trialDatum.ypos-trialDatum.ypos(1);
    if ~isempty(find(strcmp(fieldnames(bouts),'turn_power'), 1))
        trialDatum.turn_power = bouts.turn_power(:,bouts.bout_start>=trialDatum.trial_start & bouts.bout_start<= trialDatum.trial_end);
    end
    trialData{i} = trialDatum;
end
disp('Succeeded');