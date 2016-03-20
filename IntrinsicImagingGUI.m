function handles=IntrinsicImagingGUI(parent,handles)

handles.imaqGUI.hFig = figure('Toolbar','none',...
        'UserData',parent,...
       'Menubar', 'none',...
       'NumberTitle','Off',...
       'Name','Intrinsic Imaging',...
       'Position', [200   150   800   800]);
   
    % Create a figure window. This example turns off the default
    % toolbar and menubar in the figure.
    imaqreset
    adaptor_info=imaqhwinfo;
    [adaptorNum,selected]=listdlg('ListString',adaptor_info.InstalledAdaptors,...
        'SelectionMode','single',...
        'InitialValue',5,...
        'Name','Select a video adaptor:');
    if ~selected
        close_Callback(handles.imaqGUI.hFig,handles);
        return;
    end
    formats=imaqhwinfo(adaptor_info.InstalledAdaptors{adaptorNum});
    [formatNum,selected]=listdlg('ListString',formats.DeviceInfo(1).SupportedFormats,...
        'SelectionMode','single',...
        'InitialValue',3,...
        'Name','Select a video format:');
    if ~selected
        close_Callback(handles.imaqGUI.hFig,handles);
        return;
    end
    handles.vid = videoinput(adaptor_info.InstalledAdaptors{adaptorNum},...
        formats.DeviceIDs{1},...
        formats.DeviceInfo(1).SupportedFormats{formatNum});
    pause(0.5);
    handles.vid_src = getselectedsource(handles.vid);
    handles.imagingInfo.exposure=.1;
    handles.imagingInfo.satcutoff=0.90;
    handles.imagingInfo.stimf=1/8;
    handles.imagingInfo.fgi=4;
    handles.imagingInfo.frames=3000;
    handles.imagingInfo.numberPresentations=10;
    handles.imagingInfo.greenImage=[];
    handles.imagingInfo.saveVideo=4;
    handles.imagingInfo.framesize=fliplr(get(handles.vid,'VideoResolution'));

    set(handles.vid_src,'Exposure',handles.imagingInfo.exposure)
    set(handles.vid_src,'HighSensitivity','Off')
    pause(0.5)
    imaqmem(5e10);
    pause(0.5)
    set(handles.vid,'FramesPerTrigger',handles.imagingInfo.frames);
    set(handles.vid,'FrameGrabInterval',handles.imagingInfo.fgi)
   

% Set up the push buttons
uicontrol(handles.imaqGUI.hFig,'String', 'Start Preview',...
    'UserData',parent,...
    'Callback', @startPreview_Callback,...
    'Interruptible','off',...
    'Units','normalized',...
    'Position',[0.01 0.09 0.15 .07]);
uicontrol(handles.imaqGUI.hFig,'String', 'Stop Preview',...
    'UserData',parent,...
    'Callback', @stopPreview_Callback,...
    'Interruptible','off',...
    'Units','normalized',...
    'Position',[0.01 0.01 .15 .07]);
uicontrol(handles.imaqGUI.hFig,'String', 'Close',...
    'UserData',parent,...
    'Callback', @close_Callback,...
    'Interruptible','off',...
    'Units','normalized',...
    'Position',[0.84 0.01 .15 .07]);

uicontrol(handles.imaqGUI.hFig,'String', 'Take Green Snap',...
    'UserData',parent,...
    'Callback',@greenSnap_Callback,...
    'Interruptible','off',...
    'Units','normalized',...
    'Position',[.17 0.09 0.15 .07]);

uicontrol(handles.imaqGUI.hFig,'style','text','String','Exposure (secs)', ...
    'Units','normalized',...,
    'BackgroundColor',[.8 .8 .8],...
    'Position',[0.33 .15 .15 .03]);
handles.imaqGUI.hexp=uicontrol(handles.imaqGUI.hFig,'style', 'edit',...
    'string',num2str(handles.imagingInfo.exposure),...
    'Units','normalized',...
    'Position',[0.33 0.09 .15 .07]);
set(handles.imaqGUI.hexp,...
    'UserData',parent,...
    'Callback', @exposure_Callback,...
    'Interruptible','off');

uicontrol(handles.imaqGUI.hFig,'style','text','String','Cutoff', ...
    'Units','normalized',...,
    'BackgroundColor',[.8 .8 .8],...
    'Position',[0.50 .15 .15 .03]);
handles.imaqGUI.hsat=uicontrol(handles.imaqGUI.hFig,'style', 'edit',...
    'string',num2str(handles.imagingInfo.satcutoff),...
    'Units','normalized',...
    'Position',[0.50 0.09 .15 .07]);
set(handles.imaqGUI.hsat,...
    'UserData',parent,...
    'Callback', @satCutoff_Callback,...
    'Interruptible','off');

uicontrol(handles.imaqGUI.hFig,'style','text','String','Frame Grab Int.', ...
    'Units','normalized',...,
    'BackgroundColor',[.8 .8 .8],...
    'Position',[0.67 .15 .15 .03]);
handles.imaqGUI.hfgi=uicontrol(handles.imaqGUI.hFig,'style', 'edit',...
    'string',int2str(handles.imagingInfo.fgi),...
    'Units','normalized',...
    'Position',[0.67 0.09 .15 .07]);
set(handles.imaqGUI.hfgi,...
    'UserData',parent,...
    'Callback', @grabInterval_Callback,...
    'Interruptible','off');

uicontrol(handles.imaqGUI.hFig,'style','text','String','Video Import Method:', ...
    'Units','normalized',...,
    'BackgroundColor',[.8 .8 .8],...
    'Position',[0.33 .01 .15 .03]);
handles.imaqGUI.hsv=uicontrol(handles.imaqGUI.hFig,'style', 'popupmenu',...
    'string',{'Save Full Trials','Save Iterations','Save Trial Means','Save to Disk'},...
    'Units','normalized',...
    'Value',handles.imagingInfo.saveVideo,...
    'Position',[0.5 0.01 .15 .07]);
set(handles.imaqGUI.hsv,...
    'UserData',parent,...
    'Callback', @saveVideo_Callback,...
    'Interruptible','off');

% Create the text label for the timestamp
handles.imaqGUI.hTextLabel = uicontrol(handles.imaqGUI.hFig,'style','text','String','Timestamp', ...
    'Units','normalized',...
    'BackgroundColor',[.8 .8 .8],...
    'Position',[0.84 0.09 .15 .04]);

% Create the image object in which you want to
% display the video preview data.
imWidth = 696;
imHeight = 520;
nBands = get(handles.vid, 'NumberOfBands');
handles.imaqGUI.hImage = image( zeros(imHeight, imWidth, nBands) );

cmap=colormap('jet');
cmap(57:end,2:3)=0;
colormap(cmap);




% Specify the size of the axes that contains the image object
% so that it displays the image at the right resolution and
% centers it in the figure window.
figSize = get(handles.imaqGUI.hFig,'Position');
figWidth = figSize(3);
figHeight = figSize(4);
set(gca,'unit','pixels',...
        'position',[ ((figWidth - imWidth)/2)... 
                     ((figHeight - imHeight)/2+70)...
                       imWidth imHeight ]);

% Set up the update preview window function.
setappdata(handles.imaqGUI.hImage,'UpdatePreviewWindowFcn',@intrinsic_imaging_preview_fcn);

% Make handle to text label available to update function.
setappdata(handles.imaqGUI.hImage,'HandleToTimestampLabel',handles.imaqGUI.hTextLabel);

% Make handle to sat cutoff text label available to update function.
setappdata(handles.imaqGUI.hImage,'HandleToSatCutoff',handles.imaqGUI.hsat);

% preview(handles.vid);
preview(handles.vid,handles.imaqGUI.hImage);

handles.imagingInfo.previewopened=1;
guidata(handles.output,handles);
end

function close_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    if isfield(handles,'vid') && isa(handles.vid,'videoinput')
    stoppreview(handles.vid);
    pause(0.5);
    stop(handles.vid);
    pause(0.5);
    flushdata(handles.vid);
    pause(0.5);
    delete(handles.vid);
    end
    handles.imagingInfo.previewopened=0;
    set(handles.intrinsicImaging,'Value',0,'Enable','off');
    set(handles.startCamera,'Enable','on');
    rmfield(handles,'imaqGUI');
    guidata(get(hObject,'UserData'),handles);
    close(gcf);
    pause(0.5);
end

function startPreview_Callback(hObject,eventdata)
    pause(0.5);
    handles=guidata(get(hObject,'UserData'));
    preview(handles.vid);
    pause(0.5);
end

function stopPreview_Callback(hObject,eventdata)
    pause(0.5);
    handles=guidata(get(hObject,'UserData'));
    stoppreview(handles.vid);
    pause(0.5);
end

function greenSnap_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    stoppreview(handles.vid);
    handles.vid.LoggingMode='memory';
    set(handles.vid,'FramesPerTrigger',1);
    start(handles.vid);
    while isrunning(handles.vid)
        pause(0.1);
    end
    handles.imagingInfo.greenImage=squeeze(getdata(handles.vid));
    flushdata(handles.vid);
    %set(handles.vid,'FramesPerTrigger',handles.imagingInfo.frames);
    guidata(get(hObject,'UserData'),handles);
end

function numberPresentations_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    handles.imagingInfo.numberPresentations=str2num(get(handles.imaqGUI.hnp,'string'));
    guidata(get(hObject,'UserData'),handles);
end

function exposure_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    handles.imagingInfo.exposure=str2num(get(handles.imaqGUI.hexp,'string'));
    set(handles.vid_src,'exposure',handles.imagingInfo.exposure);
    guidata(get(hObject,'UserData'),handles);
end

function satCutoff_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    handles.imagingInfo.satcutoff=str2num(get(handles.imaqGUI.hsat,'string'));
    guidata(get(hObject,'UserData'),handles);
end

function movieFrames_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    handles.imagingInfo.frames=str2num(get(handles.imaqGUI.hnf,'string'));
    set(handles.vid,'FramesPerTrigger',handles.imagingInfo.frames);
    guidata(get(hObject,'UserData'),handles);
end

function grabInterval_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    handles.imagingInfo.fgi=str2num(get(handles.imaqGUI.hfgi,'string'));
    set(handles.vid,'FrameGrabInterval',handles.imagingInfo.fgi);
    guidata(get(hObject,'UserData'),handles);
end

function saveVideo_Callback(hObject,eventdata)
    handles=guidata(get(hObject,'UserData'));
    handles.imagingInfo.saveVideo=get(handles.imaqGUI.hsv,'value');
    guidata(get(hObject,'UserData'),handles);
end