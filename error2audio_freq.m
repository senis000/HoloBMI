function freq = error2audio_freq(E2mE1_error, E1_error, E2_error, fb_cal)
% Parameters from fb_cal:
% fb_cal.obj_max         = ...
%     prctile(obj_obs, fb_settings.obj_max_perctile); 
% 
% %for fb, obj is floored to this value
% fb_cal.obj_min         = ...
%     0; %for fb, cursor is floor to this value
% fb_cal.obj_range       = ...
%     fb_cal.obj_max - fb_cal.obj_min; 
% % % freq = a*exp(b*(cursor_trunc-cursor_min))
% fb_cal.a                  = ...
%     fb_cal.settings.freq_min; 
% fb_cal.b                  = ...
%     (log(fb_cal.settings.freq_max) - log(fb_cal.a))/fb_cal.obj_range; 

%z-score each error: 
% fb_cal.error_E2mE1_mu       = error_E2mE1_mu; 
% fb_cal.error_E2mE1_sigma    = error_E2mE1_sigma; 
% 
% fb_cal.error_E1_mu          = error_E1_mu; 
% fb_cal.error_E1_sigma       = error_E1_sigma; 
% 
% fb_cal.error_E2_mu          = error_E2_mu; 
% fb_cal.error_E2_sigma       = error_E2_sigma; 

error_E2mE1_z   = zscore_helper(E2mE1_error, fb_cal.error_E2mE1_mu, fb_cal.error_E2mE1_sigma);
error_E1_z      = zscore_helper(E1_error, fb_cal.error_E1_mu, fb_cal.error_E1_sigma);
error_E2_z      = zscore_helper(E2_error, fb_cal.error_E2_mu, fb_cal.error_E2_sigma);

obj = ...
    fb_cal.settings.lambda_E2mE1*error_E2mE1_z + ...
    fb_cal.settings.lambda_E1*error_E1_z + ...
    fb_cal.settings.lambda_E2*error_E2_z;

%Map objective to feedback: 


if fb_cal.settings.target_low_freq == 0
    %This means objective decreasing -> freq increasing
    %Flip the objective
    obj      = -obj;
    obj_min  = -fb_cal.obj_max;
    obj_max  = -fb_cal.obj_min;    
else
    obj_min  = fb_cal.obj_min;
    obj_max  = fb_cal.obj_max;    
end
obj_trunc      = max(obj, obj_min); 
obj_trunc      = min(obj_trunc, obj_max); 
% cursor_trunc    = max(cursor, obj_min); 
% cursor_trunc    = min(cursor_trunc, cursor_max); 

freq = fb_cal.a*exp(fb_cal.b*(obj_trunc-obj_min));
freq = double(freq); 

end

function zdata = zscore_helper(data, mu, sigma)
    zdata = (data-mu)/sigma;
end