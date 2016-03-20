function analyzePhase(datedir, experimentNumber) 
filebase=fullfile('StimGen_Results',datedir,strcat('Experiment_',int2str(experimentNumber)));
load(fullfile(filebase,strcat('imagingInfo_',datedir,'_',int2str(experimentNumber))),...
    'stimf','exposure','fgi','numberTrials','numberIters','numberConditions','imagingFrames','experimentalNotes','greenImage','saveVideo');
for trial=1:numberTrials
    for c=1:totalConditions
        switch saveVideo
            case 1
                load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                    strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))));
                [X,Y,T]=size(ISData);
                cyclestarts=startframes;
                cyclestarts(cyclestarts==0)=[];
                cyclestarts(cyclestarts==cyclestarts(end))=[];
                Fc=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
                if cyclestarts(end)+ceil(Fc)>T
                    cyclestarts(end)=[];
                    Fc=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
                end
                tcData=zeros(X,Y,ceil(Fc));
                for t=0:ceil(Fc)-1
                    tcData(:,:,t+1)=mean(ISdata(:,:,cyclestarts+t),3);
                end
            case 2
                load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                    strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'_iter_',int2str(1))));
                [X,Y,Fc]=size(ISData);
                tcData=double(ISData);
                for iter=2:numberIters
                    load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                        strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'_iter_',int2str(iter))));
                    tcData=tcData+double(ISData);
                end
            case 3
                load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                    strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))));
                load(fullfile('PSC_Data',datedir,strcat('Experiment',int2str(experimentNumber)),strcat('Condition',int2str(c)),...
                    strcat('Data_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))));
                tcData=ISdata;
                [X,Y,Fc]=size(tcData);
        end
    Fs=1/exposure/fgi; % Frames per second
    fftL=Fc;
    NFFT=round(3*Fc);%2^(nextpow2(fftL)); % 2*Fc; % 
    f=Fs/2*linspace(0,1,NFFT/2+1);
    targetf=Fs/Fc;
    fn=find(abs(f-targetf)==min(abs(f-targetf)));
    deltaF=bsxfun(@minus,tcData,mean(tcData,3));
    deltaF=smooth3(deltaF,'gaussian',[5,5,1],3);
    deltaF=smooth3(deltaF,'gaussian',[1,1,17],7);   
    avgwind=Fc*1.5+1;
    FTfn=zeros(X,Y);

    for x=1:X
        FT=fft(deltaF(x,:,:),NFFT,3)/fftL;
        FTfn(x,:)=mean(FT(:,:,fn-1:fn+1),3);
    end
    center=mean(mean(deltaF(floor(X/2-2):ceil(X/2+2),floor(Y/2-2):ceil(Y/2+2),:),1),2);
    singFT=fft(squeeze(mean(mean(center,1),2)),NFFT,1)/fftL;
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
                phasepowergreen(x,y,:)=greenmap(x,y)*cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
            end
        end
    end

    save(fullfile('PSC_Data',datedir,strcat('Experiment',int2str(experimentNumber)),strcat('Condition',int2str(c)),...
        strcat('Phase_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))),...
        'datedir','FOIphase','FOIpower','singFT' ,'-v7.3')


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

    saveas(hAnalysis,fullfile('PSC_Data',datedir,strcat('Experiment',int2str(experimentNumber)),strcat('Condition',int2str(c)),...
        strcat('PhaseFig_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'.fig')));
    saveas(hAnalysis,fullfile('PSC_Data',datedir,strcat('Experiment',int2str(experimentNumber)),strcat('Condition',int2str(c)),...
        strcat('PhaseFig_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'.jpg')));    
    
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
    
    if exist('mov')
        movie2avi(mov,fullfile('PSC_Data',datedir,strcat('Experiment',int2str(experimentNumber)),strcat('Condition',int2str(c)),...
            strcat('PhaseMov',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'.avi')),'fps',10);
    end

    end
end