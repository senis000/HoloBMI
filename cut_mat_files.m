function data=cut_mat_files(data, length_cut)
% function to cut the file short but keeping the same structure
    data.cursor(length_cut:end) = NaN;
    data.bmiAct(:, length_cut:end) = NaN;
    data.baseVector(:, length_cut:end) = NaN;
    data.selfHits(length_cut:end) = NaN;
    data.holoHits(length_cut:end) = NaN;
    data.selfVTA(length_cut:end) = NaN;
    data.holoVTA(length_cut:end)= NaN;
    data.holoDelivery(length_cut:end) = NaN;
    data.trialStart(length_cut:end) = NaN;
    data.timeVector(length_cut:end) = NaN;

    data.selTargetCounter = sum(data.selfHits);
    data.holoTargerCounter = sum(data.holoHits);
    data.selTargetVTACounter = sum(data.selfVTA);
    data.holoTargerVTACounter = sum(data.holoVTA);
    data.trialCounter = sum(data.trialStart);
    data.frame = length_cut;
end