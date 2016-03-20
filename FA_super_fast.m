for trial=1:numberTrials
    for c=1:totalConditions

    Fs=1/exposure/fgi; % Frames per second
        [X,Y,T]=size(ISdata);

    if ~isempty('startframes')
        cyclestarts=startframes;
        cyclestarts(cyclestarts==0)=[];
        cyclestarts(cyclestarts==cyclestarts(end))=[];
        Fc=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
        if cyclestarts(end)+ceil(Fc)>T
            cyclestarts(end)=[];
            Fc=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
        end
        analyzestimf=0;
    else
        Fc=Fs/stimf; % frames per cycle
        analyzestimf=1;
    end
    
    fftL=Fc;
    NFFT=round(3*Fc);%2^(nextpow2(fftL)); % 2*Fc; % 
    f=Fs/2*linspace(0,1,NFFT/2+1);

    
    targetf=Fs/Fc;

    fn=find(abs(f-targetf)==min(abs(f-targetf)));
%    cyclestarts(end-1:end)=[];
%     cyclestarts(1:3)=[];
%     temp=squeeze(mean(mean(ISdata,1),2));
%     noisysignal=find(abs(temp-mean(temp))>std(temp)*3);
%     noisycycles=zeros(size(cyclestarts));
%     for n=1:length(noisysignal)
%         cyclenum=find(noisysignal(n)>cyclestarts,1,'first');
%         if ~isempty(cyclenum) && noisysignal(n) <= cyclestarts(end)+ ceil(Fc)
%             noisycycles(cyclenum)=1;
%         end
%     end
%     cyclestarts(find(noisycycles))=[];
    cycleF=zeros(X,Y,ceil(Fc));
    for t=0:ceil(Fc)-1
        cycleF(:,:,t+1)=mean(ISdata(:,:,cyclestarts+t),3);
    end
    cycleF=repmat(cycleF,[1,1,3]);
%         deltaF=cycleF;

    deltaF=bsxfun(@minus,cycleF,mean(cycleF,3));
%     deltaF=bsxfun(@rdivide,deltaF,mean(cycleF,3));
%%
%     temp=(mean(abs(deltaF),3));
%     temp2=(mean(cycleF,3));
%     p=polyfit(temp2(:),temp(:),1)
%     deltaF=bsxfun(@rdivide,deltaF,p(1)*temp2-p(2));
% temp2=reshape(cycleF,X*Y,1,[]);
% 
%     for t=1-ceil(Fc):2*ceil(Fc)
%         temp=deltaF(:,:,t+ceil(Fc));
%         deltaF(:,:,t+ceil(Fc))=temp-mean(temp(:));
%     end
%     for t=1:ceil(Fc)
%         temp=deltaF(:,:,t);
%         deltaF(:,:,t)=(temp./prctile(abs(temp(:)),1));
%     end
%     deltaF=bsxfun(@minus,deltaF,mean(deltaF,3));
    deltaF=smooth3(deltaF,'gaussian',[5,5,1],3);
    deltaF=smooth3(deltaF,'gaussian',[1,1,17],7);
%     for t=1:ceil(Fc)
%         temp=deltaF(:,:,t);
%         deltaF(:,:,t)=(temp-mean(temp(:)));
%     end
% deltaF=bsxfun(@minus,deltaF,mean(mean(deltaF,2),1));
   
    avgwind=Fc*1.5+1;
    FTfn=zeros(X,Y);

    for x=1:X
        FT=fft(deltaF(x,:,:),NFFT,3)/fftL;
        FTfn(x,:)=mean(FT(:,:,fn-1:fn+1),3);
    end
    center=mean(mean(deltaF(floor(X/2-2):ceil(X/2+2),floor(Y/2-2):ceil(Y/2+2),:),1),2);
    singFT=fft(squeeze(mean(mean(center,1),2)),NFFT,1)/fftL;
    'fouriered'

    hasgreen=~isempty(greenImage);
    if hasgreen
        greenmap=double(greenImage)/max(double(greenImage(:)));
    end

    hAnalysis=figure;
    cmap=colormap(hsv);

FTfnsmoothed=shiftdim(smooth3(shiftdim(FTfn,-1),'gaussian',[1,7,7],5));
FOIphase=mod(angle(FTfnsmoothed)+1.5,2*pi)-pi;
    FOIpower=(abs(FTfnsmoothed));%power actually means amplitude here. real power should be square of amplitude.
maxpower=max(FOIpower(:));1e-3;
minpower=0;
    FOIphasemap=ceil((FOIphase+pi)/2/pi*64);
    FOIpowermap=max(min((FOIpower-minpower)/(maxpower-minpower),1),0);%prctile(FOIpower(:),99.9),1);
    for x=1:X
        for y=1:Y
            phasepower(x,y,:)=cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
            if hasgreen
%                 phasegreen(x,y,:)=cmap(FOIphasemap(x,y),:)/2+greenmap(x,y)/2;
                phasepowergreen(x,y,:)=greenmap(x,y)*cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
            end
        end
    end

    subplot(2,3,1)
    loglog(f,abs(squeeze(singFT(1:NFFT/2+1))).^2,'k')
    hold on
    scatter(f(fn),abs(singFT(fn)).^2,'og');
    hold off

    subplot(2,3,4)
    imagesc(repmat(FOIpowermap,[1,1,3]))
    subplot(2,3,5)
    imagesc(FOIphasemap);
    subplot(2,3,6)
    imagesc(phasepower)
    if hasgreen
        subplot(2,3,2)
        imagesc(repmat(greenmap,[1,1,3]));
        subplot(2,3,3)
        imagesc(phasepowergreen*2/3+repmat(greenmap,[1,1,3])/3)
    end

    hmov=figure;
    colormap gray;
    if exist('mov');
        clear mov
    end

    greenmax=double(max(greenImage(:)));
    overlay=zeros([size(greenImage),ceil(Fc)]);
    movieF=smooth3(deltaF,'gaussian',[3,3,3],1);
    clim=[min(movieF(:)),max(movieF(:))]*.5;
    for t=1:ceil(Fc)-1
        avgimg = deltaF(:,:,t);%mean(movieF(:,:,t+ceil(Fc)-1:t+ceil(Fc)+1),3);
        imagesc(avgimg)
        caxis(clim);
        drawnow
        pause(0.01)
        mov(t)=getframe;
        overlay(:,:,t)=max(min(avgimg,clim(2)),clim(1))+5*double(greenImage)/greenmax;
    end

%     hmov2=figure;
%     colormap gray
%     clim=[min(overlay(:)),max(overlay(:))];
%     for t=1:ceil(Fc)
%         imagesc(overlay(:,:,t))
%         caxis(clim);
%         drawnow
%         pause(.015);
%         overlaymov(t)=getframe;
%     end
toc

%{
figure
imagesc(greenimage)
BW=double(roipoly);
BW(~BW)=NaN;

    avgwind=Fc*1.5+1;
    FTfn=zeros(X,Y);

    for x=1:X
        FT=fft(deltaF(x,:,:),NFFT,3)/fftL;
        FTfn(x,:)=mean(FT(:,:,fn-1:fn+1),3);
    end
    center=mean(mean(deltaF(floor(X/2-2):ceil(X/2+2),floor(Y/2-2):ceil(Y/2+2),:),1),2);
    singFT=fft(squeeze(mean(mean(center,1),2)),NFFT,1)/fftL;
    'fouriered'

    hasgreen=~isempty(greenImage);
    if hasgreen
        greenmap=double(greenImage)/max(double(greenImage(:)));
    end

    hAnalysis=figure('Position',[60,425,1201,567]);
    cmap=colormap(hsv);

FTfnsmoothed=shiftdim(smooth3(shiftdim(FTfn,-1),'gaussian',[1,7,7],5));
FOIphase=mod(angle(FTfnsmoothed)+1.5,2*pi)-pi;
    FOIpower=(abs(FTfnsmoothed));%power actually means amplitude here. real power should be square of amplitude.
maxpower=max(FOIpower(BW(:)==1));1e-3;
minpower=0;
    FOIphasemap=ceil((FOIphase+pi)/2/pi*64);
    FOIpowermap=max(min((FOIpower-minpower)/(maxpower-minpower),1),0);%prctile(FOIpower(:),99.9),1);
    for x=1:X
        for y=1:Y
            phasepower(x,y,:)=cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
            if hasgreen
%                 phasegreen(x,y,:)=cmap(FOIphasemap(x,y),:)/2+greenmap(x,y)/2;
                phasepowergreen(x,y,:)=greenmap(x,y)*cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
            end
        end
    end

    subplot(2,3,1)
    loglog(f,abs(squeeze(singFT(1:NFFT/2+1))).^2,'k')
    hold on
    scatter(f(fn),abs(singFT(fn)).^2,'og');
    hold off

    subplot(2,3,4)
    imagesc(repmat(FOIpowermap.*BW,[1,1,3]))
    axis tight
    subplot(2,3,5)
    imagesc(FOIphasemap.*BW);
        axis tight

    subplot(2,3,6)
    imagesc(phasepower.*repmat(BW,[1,1,3]))
        axis tight

    if hasgreen
        subplot(2,3,2)
        imagesc(repmat(greenmap.*BW,[1,1,3]));
            axis tight

        subplot(2,3,3)
        imagesc((phasepowergreen+repmat(greenmap,[1,1,3])/3).*repmat(BW,[1,1,3]))
        axis tight

    end

    hmov=figure;
    colormap gray;
    if exist('mov');
        clear mov
    end

    greenmax=double(max(greenImage(:)));
    overlay=zeros([size(greenImage),ceil(Fc)]);
    movieF=smooth3(deltaF,'gaussian',[3,3,3],1);
    clim=[min(movieF(:)),max(movieF(:))]*.9;
    for t=1:ceil(Fc)
        avgimg = deltaF(:,:,t);%mean(movieF(:,:,t+ceil(Fc)-1:t+ceil(Fc)+1),3);
        imagesc(avgimg.*BW)
            axis tight

        caxis(clim);
        drawnow
        pause(0.01)
        mov(t)=getframe;
        overlay(:,:,t)=max(min(avgimg,clim(2)),clim(1))+5*double(greenImage)/greenmax;
    end
%}