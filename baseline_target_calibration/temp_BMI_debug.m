
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 

%%
sum(~isnan(data.cursor))
sum(~isnan(data.bmiAct(1,:)))
%%
valid_idxs = find(~isnan(data.cursor));
cursor_valid = data.cursor(valid_idxs); 
n_valid = zscore(data.bmiAct(:, valid_idxs),0,2); 
% dffz_valid = data.bmidffz(:, valid_idxs); 
%%
cursor_valid(1)

%%
h = plot_E_activity(n_valid', bData.E_id, E_color);
xlabel('frame'); 
ylabel('dff_z');    
title('zscore dff'); 

%%
h = figure; hold on;
offset = 3; 
for i =1:8 
    plot(valid_idxs, n_valid(i,:)-offset*i); 
end
hold on;  vline(vectorHolo(1), 'k'); hline(bData.T1); hold off

%vline(find(data.holoVTA));

% %%
% h = figure; hold on;
% offset = 3; 
% for i =1:8 
%     plot(valid_idxs, dffz_valid(i,:)-offset*i); 
% end
% hold on; vline(find(data.holoVTA)); vline(vectorHolo(1), 'k'); hline(bData.T1); hold off


%%
cursor_est = n_valid'*bData.decoder;

h = figure;
plot(cursor_est); 
hold on; vline(vectorHolo(1)); hline(bData.T1); hold off

title('estimated cursor'); 

%%
h = figure; 
plot(valid_idxs, cursor_valid); 
auxselfhits = find(data.selfHits);
hline(bData.T1); 
% hold on; vline(find(data.holoVTA)); vline(vectorHolo(1:20), 'k'); hline(bData.T1); hold off

title('real cursor');

%%
sum(cursor_valid >= bData.T1)


%%
corr(cursor_est(:), cursor_valid(:))

%%
h = figure;
hist(cursor_valid, 50); 

%%
h = figure;
hist(cursor_valid(1:500), 100); 

%%
h = figure;
hist(cursor_valid((end-500):end), 50); 


% %%
% [dff_z, cursor, target_hit, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
%     dff2cursor_target(dff, bData)