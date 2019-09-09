function roi_data = init_roi_data(im_bg, num_chan, chan_labels)
%im_mask: num_rows X num_cols x num_chan
if(numel(size(im_bg)))
    roi_data.im_bg_one_chan = im_bg; 
    im_bg = repmat(im_bg, [1 1 3]); 
end
roi_data.im_bg          = im_bg; 
roi_data.im_roi         = im_bg; 
roi_data.im_roi_rg      = im_bg; 
roi_data.num_rows       = size(im_bg,1);
roi_data.num_cols       = size(im_bg,2);
roi_data.num_rois       = 0;  
roi_data.roi_bin_cell   = {};
roi_data.x = [];
roi_data.y = [];
roi_data.r = [];
roi_data.roi_mask       = zeros(roi_data.num_rows , roi_data.num_cols); 
roi_data.roi_mask_bin   = zeros(roi_data.num_rows , roi_data.num_cols); 
roi_data.chan_logical   = []; %num_chan x num_roi

roi_data.chan = repmat(struct(...
    'label', '', ...
    'num_rois', '', ...
    'idxs', '', ...
    'im_roi',       zeros(roi_data.num_rows, roi_data.num_cols), ...
    'roi_mask',     zeros(roi_data.num_rows, roi_data.num_cols), ...
    'roi_mask_bin', zeros(roi_data.num_rows, roi_data.num_cols)), [num_chan 1]); 

for chan_i = 1:num_chan
    roi_data.chan(chan_i).label = chan_labels{chan_i}; 
end