
%%
addpath(genpath('/Users/vivekathalye/Dropbox/Code/export_fig/export_fig_11.19.15'))


%%
num_revolutions = 10
theta = repmat(linspace(0, 2*pi, 1000), [1 num_revolutions]); 
size(theta)
r = linspace(0,1,length(theta)); 
x = r.*cos(theta); 
y = r.*sin(theta); 

h = figure;
plot(x,y,'r'); 
set(gca, 'xtick', []); 
set(gca, 'xticklabel', []); 
set(gca, 'ytick', []); 
set(gca, 'yticklabel', []); 
axis square

save_dir = '/Users/vivekathalye/Desktop'
export_fig(h, fullfile(save_dir, 'spiral.eps')); 