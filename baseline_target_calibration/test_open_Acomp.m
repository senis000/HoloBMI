% test_open_Acomp.m

Acomp_file = fullfile('/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline', 'redcomp.mat'); 
exist(Acomp_file); 

d = load(Acomp_file); 

%%

idx = 1;
im = reshape(d.AComp(:,1), 256, 256); 
im = im.'; 

h = figure;
imagesc(im); 
colorbar; 
axis square

%%
