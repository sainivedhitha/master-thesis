function varargout = SEMGeditor(varargin)
% SEMGEDITOR M-file for SEMGeditor.fig
%      SEMGEDITOR, by itself, creates a new SEMGEDITOR or raises the existing
%      singleton*.
%
%      H = SEMGEDITOR returns the handle to a new SEMGEDITOR or the handle to
%      the existing singleton*.
%
%      SEMGEDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEMGEDITOR.M with the given input arguments.
%
%      SEMGEDITOR('Property','Value',...) creates a new SEMGEDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SEMGeditor_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SEMGeditor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SEMGeditor

% Last Modified by GUIDE v2.5 15-Jun-2007 13:22:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SEMGeditor_OpeningFcn, ...
                   'gui_OutputFcn',  @SEMGeditor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before SEMGeditor is made visible.
function SEMGeditor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SEMGeditor (see VARARGIN)

% Choose default command line output for SEMGeditor
handles.output = hObject;


handles.y = varargin{1};
MUPulses = varargin{2};
handles.fsamp = varargin{3};
handles.dMUAP = 25;
handles.MUAP = {};
handles.max_muap = -inf;
handles.min_muap = inf;
handles.MUPulses={};
axes(handles.axes2); handles.Ax1 = axis;
axes(handles.axes3); handles.Ax2 = axis;
handles.plotInterval = length(MUPulses)*(handles.dMUAP+21); %  length(MUPulses)*(2*handles.dMUAP+21);
handles.startPlotInterval = 1;
handles.IntervalIncrement = 500;

% Update handles structure
guidata(hObject, handles);

for k1=1:length(MUPulses)      
    handles.MUPulses{k1} = MUPulses{k1}(find(MUPulses{k1}<length(handles.y)));
end

for k1=1:length(handles.MUPulses)   
    tmp=extractMUAPocurrences(handles.y,handles.MUPulses{k1},handles.dMUAP);
    for k2=1:length(tmp)
        handles.MUAP{k1} = median(tmp{k2});     
        [tmp2 handles.MUAPcorrection(k1)] = max(handles.MUAP{k1});
        handles.max_muap = max([handles.max_muap handles.MUAP{k1}]);
        handles.min_muap = min([handles.min_muap handles.MUAP{k1}]);
    end   
end     
handles.MUAPcorrection = handles.MUAPcorrection - handles.dMUAP;
[handles.sumMUAPtrains, tmp] = extractMUAPTrains(handles.y,handles.fsamp,handles.MUPulses,handles.dMUAP,0);

handles.max_y = max(handles.y);
handles.min_y = min(handles.y);
handles.p2p_y = handles.max_y - handles.min_y;

% plot signals & MUAP templates
plotAxes2(hObject, eventdata, handles);

axes(handles.axes2); 
handles.Ax1 = axis;

plotAxes3(hObject, eventdata, handles);


% UIWAIT makes SEMGeditor wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = SEMGeditor_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.startPlotInterval = handles.startPlotInterval + handles.IntervalIncrement;
if handles.startPlotInterval > length(handles.y) - handles.plotInterval;
    handles.startPlotInterval = length(handles.y) - handles.plotInterval - 1;
end
% Update handles structure
guidata(hObject, handles);

% plot signals & MUAP templates
plotAxes2(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.startPlotInterval = handles.startPlotInterval - handles.IntervalIncrement;
if handles.startPlotInterval < 1
    handles.startPlotInterval = 1;
end
% Update handles structure
guidata(hObject, handles);

% plot signals & MUAP templates
plotAxes2(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- replots the Axes2.
function plotAxes2(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2); 
hold off;
plot([handles.startPlotInterval : handles.startPlotInterval+handles.plotInterval]/handles.fsamp, ...
      -1.2*handles.min_y+handles.y(handles.startPlotInterval : handles.startPlotInterval+handles.plotInterval),'b','LineWidth',2);
hold on;
plot([handles.startPlotInterval : handles.startPlotInterval+handles.plotInterval]/handles.fsamp,...
      handles.y(handles.startPlotInterval : handles.startPlotInterval+handles.plotInterval) - ...
      handles.sumMUAPtrains(handles.startPlotInterval : handles.startPlotInterval+handles.plotInterval),'b','LineWidth',2);
xlabel('Time [s]','FontSize',14);
for k1 = 1:length(handles.MUPulses)
    tmp = handles.MUPulses{k1}(find(handles.MUPulses{k1} > handles.startPlotInterval ...
                                    & handles.MUPulses{k1} < handles.startPlotInterval+handles.plotInterval ));
    for k2 = 1:length(tmp)        
        text( ( tmp(k2) + handles.MUAPcorrection(k1) )/handles.fsamp,0.5-0.06*k1*handles.min_y,num2str(k1),'FontSize',12);        
    end   
end
axis tight; 
handles.Ax1 = axis;
axis([handles.Ax1(1:2) 0.3*handles.min_y  -1.2*handles.min_y+handles.max_y  ]);
handles.Ax1 = axis;

% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- replots the Axes3.
function plotAxes3(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)axes(handles.axes3); 
axes(handles.axes3); cla;
hold on; 

dyLine = 1.1*ones(1,length(handles.MUAP)) * (handles.max_muap - handles.min_muap);
dyLine(ceil(length(handles.MUAP)/2)+1:end) = 0;
dxLine = [ 1:ceil(length(handles.MUAP)/2) 1:(length(handles.MUAP) - ceil(length(handles.MUAP)/2))];
for k1 = 1:length(handles.MUAP)        
    plot(dxLine(k1)*(length(handles.MUAP{k1})+20) + [1:length(handles.MUAP{k1})], dyLine(k1) + handles.MUAP{k1},'b','LineWidth',2);
    text(dxLine(k1)*(length(handles.MUAP{k1})+30)+length(handles.MUAP{k1})/100, dyLine(k1) + 0.8*handles.min_muap,num2str(k1),'FontSize',14);
end
axis tight;
handles.Ax2 = axis;

%newHeight = (handles.Ax1(4)-handles.Ax1(3)) / 3.1492;
axis([handles.Ax2(1:2) ( handles.Ax1(3:4) - handles.Ax1(4) + handles.Ax2(4) )]);
handles.Ax2 = axis;
% set(handles.axes3,'YTick',[]);
% set(handles.axes3,'YTickLabel',[]);
% set(handles.axes3,'XTick',[]);
% set(handles.axes3,'XTickLabel',[]);

% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%