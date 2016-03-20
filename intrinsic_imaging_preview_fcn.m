function intrinsic_imaging_preview_fcn(obj,event,himage)
% Example update preview window function.
% handles=guidata(get(obj,'UserData'));
% Get timestamp for frame.
 tstampstr = event.Timestamp;

% Get handle to text label uicontrol.
 ht = getappdata(himage,'HandleToTimestampLabel');
 hs = getappdata(himage,'HandleToSatCutoff');

 % Set the value of the text label.
 set(ht,'String',tstampstr);
satcutoff=str2num(get(hs,'string'));
% Display image data.
class(event.Data)
max(event.Data(:))
satevent=(1+imresize(double(event.Data),[520, 696]))/256;
satevent(satevent>0.9)=0;
satoff=satevent;
satoff(satoff>=satcutoff)=0;

set(himage, 'CData', cat(3,satevent,satoff,satoff))