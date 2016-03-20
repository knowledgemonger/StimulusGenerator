function analyzePhaseAggregate(datedir, experimentNumber,post_processing,filemodifier,numT)
if nargin<3
    post_processing=1;
end
if nargin<4
    switch post_processing
        case 1
            filemodifier='dFoverf_';
        case 2
            filemodifier='dF_';
        case 3
            filemodifier='dFoverF_norm_';
        case 4
            filemodifier='dF_norm_';
        case 5
            filemodifier='dNoverF_';
    end
end
if nargin<5
    numT=0;
end
filebase=fullfile('StimGen_Results',datedir,strcat('Experiment_',int2str(experimentNumber)));
load(fullfile(filebase,strcat('imagingInfo_',datedir,'_',int2str(experimentNumber))),...
    'stimf','exposure','fgi','numberTrials','numberIters','numberConditions','imagingFrames','experimentalNotes','greenImage','saveVideo');
totalConditions=prod(numberConditions);
if numT
    numberTrials=numT;
end
for c=1:totalConditions
        aggData{c}=0;
        for trial=1:numberTrials
            switch saveVideo
                case 1
                    load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                        strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))));
                    [X,Y,T]=size(ISdata);
                    cyclestarts=startframes+1;
%                     cyclestarts(cyclestarts==0)=[];
%                     cyclestarts(cyclestarts==cyclestarts(end))=[];
                    if trial==1
                        Fc(c)=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
                    end
                    if cyclestarts(end)+ceil(Fc(c))>T
                        cyclestarts(end)=[];
                    end
                    if trial==1
                        Fc(c)=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
                    end
                    tcData=zeros(X,Y,ceil(Fc(c)));
                    for t=0:ceil(Fc(c))-1
                        tcData(:,:,t+1)=mean(double(ISdata(:,:,cyclestarts+t)),3);
                    end
                case 4
                    load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                        strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))));
                    vidreader=VideoReader(logfile);
                    ISdata=read(vidreader);
                    [X,Y,T]=size(ISdata);
                    cyclestarts=startframes+1;
%                     cyclestarts(cyclestarts==0)=[];
%                     cyclestarts(cyclestarts==cyclestarts(end))=[];
                    if trial==1
                        Fc(c)=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
                    end
                    if cyclestarts(end)+ceil(Fc(c))>T
                        cyclestarts(end)=[];
                    end
                    if trial==1
                        Fc(c)=(cyclestarts(end)-cyclestarts(1))/(length(cyclestarts)-1);
                    end
                    tcData=zeros(X,Y,ceil(Fc(c)));
                    for t=0:ceil(Fc(c))-1
                        tcData(:,:,t+1)=mean(double(ISdata(:,:,cyclestarts+t)),3);
                    end
                case 2
                    load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                        strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'_iter_',int2str(1))));
                    [X,Y,Fc(c)]=size(ISdata);
                    tcData=double(ISdata);
                    for iter=2:numberIters
                        load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                            strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial),'_iter_',int2str(iter))));
                        tcData=tcData+double(ISdata);
                    end
                case 3
                    load(fullfile(filebase,strcat('Condition_',int2str(c)),...
                        strcat('imagingData_',datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'_trial_',int2str(trial))));
                    tcData=ISdata;
                    [X,Y,Fc(c)]=size(tcData);
            end
            aggData{c}=aggData{c}+tcData/numberTrials;
        end
end

mean_t=0;
for c=1:totalConditions
    mean_t=mean_t+mean(aggData{c},3)/totalConditions;
end

Fc
for c=1:totalConditions
    Fs=1/exposure/fgi; % Frames per second
    fftL=Fc(c);
    NFFT=round(3*Fc(c));%2^(nextpow2(fftL)); % 2*Fc(c); % 
    f=Fs/2*linspace(0,1,NFFT/2+1);
    targetf=Fs/Fc(c);
    fn=find(abs(f-targetf)==min(abs(f-targetf)));
    mean_xy=mean(mean(aggData{c},2),1);
    switch post_processing
        case 1
            deltaF=bsxfun(@rdivide,aggData{c},mean_t)-1;
            deltaF=bsxfun(@minus,deltaF,mean(deltaF,3));
        case 2
            deltaF=bsxfun(@minus,aggData{c},mean_t);
        case 3
            deltaF=bsxfun(@rdivide,aggData{c},mean_t)-1;
            deltaF=bsxfun(@rdivide,deltaF,mean_xy);
        case 4
            deltaF=bsxfun(@minus,aggData{c},mean_t);
            deltaF=bsxfun(@rdivide,deltaF,mean_xy);
        case 5
            deltaF=bsxfun(@rdivide,bsxfun(@rdivide,aggData{c},mean_t),mean_xy);
            deltaF=bsxfun(@minus,deltaF,mean(deltaF,3));
    end
%     deltaF=smooth3(deltaF,'gaussian',[5,5,7],7);
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
    FOIphase=angle(FTfnsmoothed);
    FOIpower=(abs(FTfnsmoothed));%power actually means amplitude here. real power should be square of amplitude.
    maxpower=max(FOIpower(:));1e-3;
    minpower=0;
    FOIphasemap=ceil((FOIphase+pi)/2/pi*64);
    FOIpowermap=max(min((FOIpower-minpower)/(maxpower-minpower),1),0);%prctile(FOIpower(:),99.9),1);
    phasepower=zeros(X,Y,size(cmap,2));
    phasepowergreen=phasepower;
    for x=1:X
        for y=1:Y
            if isnan(FOIphasemap(x,y))
                phasepower(x,y,:)=NaN;
                if hasgreen
                    phasepowergreen(x,y,:)=NaN;
                end
            else
                phasepower(x,y,:)=cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
                if hasgreen
                    phasepowergreen(x,y,:)=greenmap(x,y)*cmap(FOIphasemap(x,y),:)*FOIpowermap(x,y);
                end
            end
        end
    end

    save(fullfile(filebase,strcat('Condition_',int2str(c)),...
        strcat('phaseResults_',filemodifier,datedir,'_',int2str(experimentNumber),'_c_',int2str(c))),...
        'datedir','FOIphase','FOIpower','singFT', 'FTfn','-v7.3')


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

    saveas(hAnalysis,fullfile(filebase,strcat('Condition_',int2str(c)),...
        strcat('phaseFig_',filemodifier,datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'.fig')));
    saveas(hAnalysis,fullfile(filebase,strcat('Condition_',int2str(c)),...
        strcat('phaseFig_',filemodifier,datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'.jpg')));    
    
    hmov=figure;
    colormap gray;


    movieF=smooth3(deltaF,'gaussian',[3,3,3],1);
    clim=[min(movieF(:)),max(movieF(:))]*.5;
    mov(1:ceil(Fc(c))) = struct('cdata', [],'colormap', []);
    for t=1:ceil(Fc(c))
        avgimg = mean(movieF(:,:,t),3);
        imagesc(avgimg)
        caxis(clim);
        drawnow
        pause(0.01)
        mov(t)=getframe(hmov);
    end
    movie2avi(mov,fullfile(filebase,strcat('Condition_',int2str(c)),...
        strcat('phaseMov_',filemodifier,datedir,'_',int2str(experimentNumber),'_c_',int2str(c),'.avi')),'fps',5);
    end