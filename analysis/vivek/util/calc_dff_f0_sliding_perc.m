function [dff, f0] = calc_dff_f0_sliding_perc(ts, f0_perc, f0_win)

% %DEUBG: 
% ts = target_f_ts.';
% f0_win = 2*60*30; 
% f0_perc = 10; 

%%
tic
disp('calculating F0 sliding percentile: ');
f0_no_prefix = sliding_perc(ts, f0_win, f0_perc);
disp('DONE F0 sliding percentile: ');
toc

%%
f0 = [f0_no_prefix(1)*ones(f0_win-1, 1); f0_no_prefix]; %pad the beginning
dff = (ts - f0)./(f0);

%%
% % DEBUG: 
% sel = 1:length(ts); %1:5000 + 4000; 
% h = figure;
% hold on;
% plot(ts(sel)); 
% plot(f0(sel), 'Color', 'k', 'LineWidth', 2); 

%%
% %DEBUG: 
% h = figure;
% plot(dff); 
% title('dff'); 
% [dff_norm,~] = range_norm(dff, 100, 1);
% [f_norm,~] = range_norm(ts, 100, 1); 
% 
% h = figure;
% hold on; 
% plot(f_norm); 
% plot(dff_norm); 
% 
% legend('f norm','dff norm'); 