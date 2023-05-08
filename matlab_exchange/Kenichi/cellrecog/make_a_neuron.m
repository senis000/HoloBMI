function im=make_a_neuron (diam,fill);

x=diam*cos(1:fill);
y=diam*sin(1:fill);
figure(4), 
    set(gcf,'Position',[300 300 90 90]);
h=plot(x,y); 
    box off;
    set(gca,'YTickLabel',[])
    set(gca,'YTick',[])
    set(gca,'XTickLabel',[])
    set(gca,'XTick',[])
    

im=getframe;

im=horzcat(im.cdata);
im=im(2:end-1,2:end-1,1);
im(find(im==0))=1;
im(find(im==255))=0;
figure(5), imagesc(im)
