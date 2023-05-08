function [roi_t] = delete_neuropil(roi_t)


imagesc(roi_t.raw), colormap bone, caxis([0 500])
mou = drawROIs_kk(1);
neur_map = extractROIs2;
imagesc(roi_t.raw), colormap bone, caxis([0 500])
mou2 = drawROIs_kk(1);
neur_map2 = extractROIs2;

dele_mask = ((neur_map + neur_map2)>0).*roi_t.exmask;
dele_inmask = ((neur_map + neur_map2)>0).*roi_t.inmask;
to_dele = unique(dele_mask)';
into_dele = unique(dele_inmask)';
to_dele(to_dele==0) = [];
into_dele(into_dele==0) = [];
[roi_t] = delete_neurons(roi_t,to_dele,'e');
[roi_t] = delete_neurons(roi_t,into_dele,'i');

end