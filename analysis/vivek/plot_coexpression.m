function [fov_im, E2_mask, E1_mask, bmi_mask, E1_overlay, E2_overlay, BMI_overlay] = ...
    plot_coexpression(mask_rg_data, tinfo, red_scale, green_scale, plot_dir)
%6.8.19
%Assumes input:
%mask_rg_data
%                   Im: [512x512 double]
%                  Img: [512x512 double]
%             holoMask: [512x512 double]
%     holoMaskRedGreen: [512x512 double]
%                  red: [2x185 double]
%             redGreen: [2x189 double]
%
%plot_dir
%Unpack:
Im                  = mask_rg_data.Im; 
Img                 = mask_rg_data.Img; 
holoMask            = mask_rg_data.holoMask; 
holoMaskRedGreen    = mask_rg_data.holoMaskRedGreen;
red                 = mask_rg_data.red; 
redGreen            = mask_rg_data.redGreen;

Im_rescale = Im-min(Im(:));
Im_rescale = red_scale*Im_rescale/max(Im_rescale(:));

%
Img_rescale = Img-min(Img(:)); 
Img_rescale = green_scale*Img_rescale/max(Img_rescale(:)); 

%
fov_im = zeros(512,512,3); 
fov_im(:,:,1) = Im_rescale; 
fov_im(:,:,2) = Img_rescale; 
h=figure; imshow(fov_im)

export_fig(h, fullfile(plot_dir, 'coexpression.png')); 
export_fig(h, fullfile(plot_dir, 'coexpression.eps')); 

%BMI Neuron mask
E2_mask     = zeros(size(holoMaskRedGreen)); 
for indn=1:length(tinfo.E2_base)
    auxmask = holoMaskRedGreen;
    auxmask(auxmask~=tinfo.E2_base(indn)) = 0;
    auxmask(auxmask~=0) = indn;
    E2_mask = auxmask + E2_mask;    
end

E1_mask     = zeros(size(holoMaskRedGreen)); 
for indn=1:length(tinfo.E1_base)
    auxmask = holoMaskRedGreen;
    auxmask(auxmask~=tinfo.E1_base(indn)) = 0;
    auxmask(auxmask~=0) = indn;
    E1_mask = auxmask + E1_mask;    
end

%E1
h = figure;
imagesc(E1_mask); 
axis square
colorbar
title('E1 mask'); 
export_fig(h, fullfile(plot_dir, 'E1_mask.png')); 
export_fig(h, fullfile(plot_dir, 'E1_mask.eps')); 

%
%E1_overlay: 
E1_overlay = fov_im;
E1_overlay(:,:,3) = E1_mask; 
h=figure; imshow(E1_overlay);
export_fig(h, fullfile(plot_dir, 'E1_overlay.png')); 
export_fig(h, fullfile(plot_dir, 'E1_overlay.eps')); 

%E2
h = figure;
imagesc(E2_mask); 
axis square
colorbar
title('E2 mask'); 
export_fig(h, fullfile(plot_dir, 'E2_mask.png')); 
export_fig(h, fullfile(plot_dir, 'E2_mask.eps')); 

%E2_overlay: 
E2_overlay = fov_im;
E2_overlay(:,:,3) = E2_mask; 
h=figure; imshow(E2_overlay);
export_fig(h, fullfile(plot_dir, 'E2_overlay.png')); 
export_fig(h, fullfile(plot_dir, 'E2_overlay.eps')); 

%BMI:
tmp = E2_mask;
tmp(tmp>0) = tmp(tmp>0)+length(tinfo.E1_base); 
bmi_mask = E1_mask + tmp; 
h = figure;
imagesc(bmi_mask); 
axis square
colorbar
title('BMI mask'); 
export_fig(h, fullfile(plot_dir, 'BMI_mask.png')); 
export_fig(h, fullfile(plot_dir, 'BMI_mask.eps')); 

%BMI_overlay: 
BMI_overlay = fov_im;
BMI_overlay(:,:,3) = bmi_mask; 
h=figure; imshow(BMI_overlay);
export_fig(h, fullfile(plot_dir, 'BMI_overlay.png')); 
export_fig(h, fullfile(plot_dir, 'BMI_overlay.eps')); 

%BMI overlay E1, E2:
fov_im = zeros(512,512,3); 
fov_im(:,:,2) = Img_rescale; 
fov_im(:,:,1) = 0.2*E2_mask>0; 
fov_im(:,:,3) = 0.2*E1_mask>0; 
h = figure;
imshow(fov_im); 
axis square
export_fig(h, fullfile(plot_dir, 'fov_E1_blue_E2_red.png')); 
export_fig(h, fullfile(plot_dir, 'fov_E1_blue_E2_red.eps')); 

close all;