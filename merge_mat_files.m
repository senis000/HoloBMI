% function to unite two bmi_online
% each case may be different 

% animal= NVI12 date = 191007
animal= 'NVI12';
date = '191007';
folder = 'C:/Users/Nuria/Documents/DATA/holoBMI/raw_data/191007/NVI12/D12/';
data = data_1;
bData = bData_1;
data.selfTargetCounter = data_1.selfTargetCounter + data_2.selfTargetCounter;
data.holoTargetCounter = data_1.holoTargetCounter + data_2.holoTargetCounter;
data.selfTargetVTACounter = data_1.selfTargetVTACounter + data_2.selfTargetVTACounter;
data.holoTargetVTACounter = data_1.holoTargetVTACounter + data_2.holoTargetVTACounter;
data.schedHoloCounter = data_1.schedHoloCounter + data_2.schedHoloCounter;
data.schedVTACounter = data_1.schedVTACounter + data_2.schedVTACounter;
data.trialCounter = data_1.trialCounter + data_2.trialCounter - 1;

add_cursor = data_2.cursor(~isnan(data_2.cursor));
data.frame = data_1.frame + length(add_cursor);

add_bmiAct = data_2.bmiAct(:,~isnan(data_2.cursor));
add_baseVector = data_2.baseVector(:,~isnan(data_2.cursor));
add_selfHits = data_2.selfHits(~isnan(data_2.cursor));
add_holoHits = data_2.holoHits(~isnan(data_2.cursor));
add_selfVTA = data_2.selfVTA(~isnan(data_2.cursor));
add_holoVTA = data_2.holoVTA(~isnan(data_2.cursor));
add_holoDelivery = data_2.holoDelivery(~isnan(data_2.cursor));
add_trialStart = data_2.trialStart(~isnan(data_2.cursor));
add_timeVector = data_2.timeVector(~isnan(data_2.cursor));
if length(add_cursor) > length(data.cursor) - (data_1.frame -1)
    data.cursor(data_1.frame:end) = add_cursor(1:(length(data.cursor) - data_1.frame-1));
    data.bmiAct(:,data_1.frame:end) = add_bmiAct(:,1:(length(data.cursor) - data_1.frame-1));
    data.baseVector(:,data_1.frame:end) = add_baseVector(:,1:(length(data.cursor) - data_1.frame-1));
    data.selfHits(data_1.frame:end) = add_selfHits(1:(length(data.cursor) - data_1.frame-1));
    data.holoHits(data_1.frame:end) = add_holoHits(1:(length(data.cursor) - data_1.frame-1));
    data.selfVTA(data_1.frame:end) = add_selfVTA(1:(length(data.cursor) - data_1.frame-1));
    data.holoVTA(data_1.frame:end) = add_holoVTA(1:(length(data.cursor) - data_1.frame-1));
    data.holoDelivery(data_1.frame:end) = add_holoDelivery(1:(length(data.cursor) - data_1.frame-1));
    data.trialStart(data_1.frame:end) = add_trialStart(1:(length(data.cursor) - data_1.frame-1));
    data.timeVector(data_1.frame:end) = add_timeVector(1:(length(data.cursor) - data_1.frame-1));
else
    data.cursor(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_cursor;
    data.bmiAct(:,data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_bmiAct;
    data.baseVector(:,data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_baseVector;
    data.selfHits(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_selfHits;
    data.holoHits(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_holoHits;
    data.selfVTA(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_selfVTA;
    data.holoVTA(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_holoVTA;
    data.holoDelivery(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_holoDelivery;
    data.trialStart(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_trialStart;
    data.timeVector(data_1.frame:(data_1.frame+length(add_cursor)-1)) = add_timeVector;
end

file = 'BMI_merged.mat'; 
save(fullfile(folder, file), 'bData', 'data')

