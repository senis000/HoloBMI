numimgs = size(imfinfo(filename),1);
for i=1:numimgs
imag(:,:,i) = imread(filename,i);
end

numneurons = max(max(holoMaskRedGreen));
unitVals = zeros(numneurons, numimgs);
backVals = zeros(1,max(max(holoMaskRedGreen)));
masktocheck = holoMaskRedGreen;
masktocheck(masktocheck>0)=nan;
masktocheck(~isnan(masktocheck)) = 1;
for i=1:numimgs
    unitVals(:, i)=obtainRoi(imag(:,:,i), strcMask);
    backVals (i) = nanmean(nanmean(imag(:,:,i)*masktocheck));
end

numstims=20;
psth=zeros(numstims, numneurons, 3*30);
for nn=1:numstims
    psth(nn,:,:) = unitVals(:,nn*150-30:nn*150+2*30-1);
end

psthmean = zeros(4,numneurons, 3*30);
psthstd = zeros(4,numneurons);
for nn=1:4
    psthmean(nn,:,:) = nanmean(psth(nn:4:end,:,:),1);
    for neur=1:numneurons
        psthstd(nn,neur) = nanmean(nanstd(squeeze(psth(nn:4:end,neur,:)),0,1));
    end
end


labelneurons = zeros(4,numneurons);
for ll=1:numneurons
    for nn=1:4
        base = psthmean(nn,ll,1:10);
        post = nanmean(psthmean(nn,ll,11:end));
        labelneurons(nn,ll)= post > nanmean(base)+2*nanstd(base);
    end
end


[lx,ly] = find(labelneurons==1);
for nn=1:length(lx)
    subplot(3,2,nn)
    plot(-10:79,squeeze(psthmean(lx(nn),ly(nn),:)))
    xlabel(squeeze(psthstd(lx(nn),ly(nn))))
end



