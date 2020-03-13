%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is developed to simulate the diffusion phenomenon governed %
% by the Allen-Cahn equation and Cahn-Hillard equation in 1D case.       % 
%                                                                        %
% Author: Zirui Mao (Post-doc in Professor Demkowicz's group)            %
% Date last modified: March., 10th, 2020                                 %
% This script is developed and can be used only for the course 'MSEN 620 %
% KINETIC PROCESS MAT SCI' instructed by Dr. Demkowicz.                  % 
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!!!!!! Read me before using !!!!!!!!!!                              %
% By using this freeware, you are agree to the following:                %
% 1. you are free to copy and redistribute the material in any format;   %
% 2. you are free to remix, transform, and build upon the material for   %
%    any purpose, even commercially;                                     %
% 3. you must provide the name of the creator and attribution parties,   %
%    a copyright notice, a license notice, a disclaimer notice, and a    % 
%    link to the material (https://github.com/maozirui/PFM.1D);%
% 4. users are entirely at their own risk using this freeware.           %
%                                                                        %
% Before use, please read the License carefully:                         %
% <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">   %
% <img alt="Creative Commons License" style="border-width:0"             %
% src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />    %
% This work is licensed under a                                          %
% <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">   %
% Creative Commons Attribution 4.0 International License</a>.            %
%                                                                        %
% %%%%%%%%% Numerical model:                                             %
%  /|\                             |                                     %
%   |                              |                                     %
%   |                              |                                     %
%   |                                                                    %
%   |                   ___________C2___________                         %
%   |                  |           |            |                        %
%   |                  |           |            |                        %
%   |                  |                        |                        %
%   |_______C1_________|<-------- w=r*L ------->|__________C1________    %
%   |                                                                    %
%   |<---------------------------- L ------------------------------->|   %
%   |__________________________________________________________________\ %
%   0                              |                                   / %
%                                  |                                     %
% %%%%%%%%% Governing equation:                                          %
%                                                                        %
%   Allen-Cahn equation:                                                 %
%   dc/dt = - M * [\Deltaf * (df/dc) - \gamma * div(c)]                  %
%   Cahn-Hillar equation:                                                %
%   dc/dt = M * div(\Deltaf * (df/dc) - \gamma * div(c))                 %
%                                                                        %
% where c = phase field (also always indicated by \phi) c in [-1, 1]     %
%       M = Mobility coefficient                                         %
%       f = local potential energy = (c^2-1)^2                           %
%       \Deltaf = local potential energy difference                      %
%       \gamma = gradient energy coefficient                             %
%       L = length of 1D domain                                          %
%       w = width of the central segment controlled by the factor f      %
%       r = ratio of the central segment to L.  Note: r belong (0,1)     %
%      C1 = one constant defining the initial c in both sides            %
%      C2 = another constant defining the initial c in central segment   %
%       iter = total number of iterations                                %
%                                                                        %
% %%%%%%%%% Boundary Conditions:                                         %
% Periodical boundary is applied to the both ending nodes, i.e.,         %
%      x(N) | x(1) x(2) .................. x(N-1) x(N) | x(1)            %
%                                                                        %
% %%%%%%%%% Finite Difference Approximation:                             %
% The second-order derivative of c in the govering equation (1) is       %
% approximated with the 2nd order accurate central Finite Difference     %
% scheme, i.e.,                                                          %
%                                                                        %
%  d^2(c)   c(i+1) - 2c(i) + c(i-1)                                      %
%  ------ = -----------------------                                      %
%   dx^2           (dx)^2                                                %
%                                                                        %
% %%%%%%%%% inputs:                                                      %
% L, r, C1, C2; M, \gamma, \Deltaf, iter;  A-C or C-H                    %
%                                                                        %  
%                                                                        %
%%%%%%%%%% How to use this freeware:                                     %
% 1. Download the documents 'Diffusion.m' and 'Diffusion.fig';           %
% 2. Open the 'Difussion.m' with Matlab.                                 %
% 3. Before running it, please ensure the two documents locate in the    %
%    same folder.                                                        %
% 4. Click 'Run' in Matlab toolstrip under 'EDITOR'.                     %
% 5. The GUI (User Interface) shall show up.                             %
% 6. Specify the values for each user-define parameters;                 %
% 7. Click the 'Run' botton in GUI when everything is ready;             %
% 8. The instant result will be plotted in the right panel;              %
% 9. Click the 'Stop' button in GUI to stop the calculation at any time; %
% 10. An animation will be saved in the document-located folder after    %
%     each time of calculation. Please rename it immediately otherwise   %
%     it will be overwritten.                                            %
% 11. Please assign credit to https://github.com/maozirui/PFM.1D for use %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varargout = Diffusion(varargin)
% DIFFUSION MATLAB code for Diffusion.fig
%      DIFFUSION, by itself, creates a new DIFFUSION or raises the existing
%      singleton*.
%
%      H = DIFFUSION returns the handle to a new DIFFUSION or the handle to
%      the existing singleton*.
%
%      DIFFUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIFFUSION.M with the given input arguments.
%
%      DIFFUSION('Property','Value',...) creates a new DIFFUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Diffusion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Diffusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Diffusion

% Last Modified by GUIDE v2.5 09-Mar-2020 22:25:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Diffusion_OpeningFcn, ...
                   'gui_OutputFcn',  @Diffusion_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Diffusion is made visible.
function Diffusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Diffusion (see VARARGIN)

javaFrame = get(gcf,'JavaFrame');
set(javaFrame,'Maximized',1);

L=1; % Length
M=1;  % collective diffusion coefficient
df=1; % local potential difference
gamma = 0.001; % gradient energy coefficient
C1=-1; % constant 1
C2=1; % constant 2
r=0.5; % ratio of central segment to the total length
N=500; % total number of nodes
x=[0:1/N:L]; % x location
c=x; % density of diffusing material
Pl=L*(1-r)/2; % left piecewise point
Pr=L-Pl; % right piecewise point

%%%% initial state %%%%%%%
c(abs(x-L/2)<=r*L/2) = C2; 
c(abs(x-L/2)>r*L/2) = C1;
c0=c;

plot(x,c,'k','linewidth',1.0);
xlabel('x');
ylabel('\phi');
set(gca,'TickDir','out');
axis([0 L min(C1,C2)-abs(C1-C2)*0.5 max(C1,C2)+abs(C1-C2)*0.5]);
set(gca,'Fontname','Times New Roman','Fontsize',16);

% Choose default command line output for Diffusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Diffusion wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Diffusion_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Get default command line output from handles structure
varargout{1} = handles.output;


function L_input_Callback(hObject, eventdata, handles)
% hObject    handle to L_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
L=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of L_input as text
%        str2double(get(hObject,'String')) returns contents of L_input as a double


% --- Executes during object creation, after setting all properties.
function L_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r_input_Callback(hObject, eventdata, handles)
% hObject    handle to r_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_input as text
%        str2double(get(hObject,'String')) returns contents of r_input as a double
r=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function r_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C1_input_Callback(hObject, eventdata, handles)
% hObject    handle to C1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
C1=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of C1_input as text
%        str2double(get(hObject,'String')) returns contents of C1_input as a double


% --- Executes during object creation, after setting all properties.
function C1_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C2_input_Callback(hObject, eventdata, handles)
% hObject    handle to C2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
C2=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of C2_input as text
%        str2double(get(hObject,'String')) returns contents of C2_input as a double


% --- Executes during object creation, after setting all properties.
function C2_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function M_input_Callback(hObject, eventdata, handles)
% hObject    handle to M_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of M_input as text
%        str2double(get(hObject,'String')) returns contents of M_input as a double

% --- Executes during object creation, after setting all properties.
function M_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma_input_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gamma=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of gamma_input as text
%        str2double(get(hObject,'String')) returns contents of gamma_input as a double


% --- Executes during object creation, after setting all properties.
function gamma_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function df_input_Callback(hObject, eventdata, handles)
% hObject    handle to df_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
df=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of df_input as text
%        str2double(get(hObject,'String')) returns contents of df_input as a double


% --- Executes during object creation, after setting all properties.
function df_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to df_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iter_input_Callback(hObject, eventdata, handles)
% hObject    handle to iter_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iter=str2double(get(hObject,'String'));
update_button_Callback(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of iter_input as text
%        str2double(get(hObject,'String')) returns contents of iter_input as a double


% --- Executes during object creation, after setting all properties.
function iter_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eqn_type = get(handles.popupmenu1, 'value');
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
L=str2double(get(handles.L_input,'String'));
r=str2double(get(handles.r_input,'String'));
C1=str2double(get(handles.C1_input,'String'));
C2=str2double(get(handles.C2_input,'String'));
M=str2double(get(handles.M_input,'String'));
gamma=str2double(get(handles.gamma_input,'String'));
df=str2double(get(handles.df_input,'String'));
iter=str2double(get(handles.iter_input,'String'));
N=500; % total number of nodes
x=[0:1/N:L]; % x location
c=x; % density of diffusing material
Pl=L*(1-r)/2; % left piecewise point
Pr=L-Pl; % right piecewise point

%%%% initial state %%%%%%%
c(abs(x-L/2)<=r*L/2) = C2; 
c(abs(x-L/2)>r*L/2) = C1;
c0=c;

plot(x,c,'k','linewidth',1.0);
xlabel('x');
ylabel('\phi');
set(gca,'TickDir','out');
axis([0 L min(C1,C2)-abs(C1-C2)*0.5 max(C1,C2)+abs(C1-C2)*0.5]);
set(gca,'Fontname','Times New Roman','Fontsize',17);


% --- Executes on button press in Run_button.
function Run_button_Callback(hObject, eventdata, handles)
% hObject    handle to Run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
L=str2double(get(handles.L_input,'String'));
r=str2double(get(handles.r_input,'String'));
C1=str2double(get(handles.C1_input,'String'));
C2=str2double(get(handles.C2_input,'String'));
M=str2double(get(handles.M_input,'String'));
gamma=str2double(get(handles.gamma_input,'String'));
df=str2double(get(handles.df_input,'String'));
iter=str2double(get(handles.iter_input,'String'));
eqn_type = get(handles.popupmenu1, 'value');
N=500; % total number of nodes
x=[0:1/N:L]; % x location
c=x; % density of diffusing material
Pl=L*(1-r)/2; % left piecewise point
Pr=L-Pl; % right piecewise point
nprint = 100; % total number of plots
%%%% initial state %%%%%%%
c(abs(x-L/2)<=r*L/2) = C2;
c(abs(x-L/2)>r*L/2) = C1;
c0=c;

dx=1/N; % grid size
if eqn_type == 1 % A-C equation
    dt=0.2*dx*dx/2/M/gamma; % critical time step dt by following the von Neumann criterion
elseif eqn_type == 2 % C-H equation
    dt = dx^2/M/(4+8*gamma/dx/dx)*0.2; % critical time 
end

nt=iter; % total time steps
nprint=max(floor(iter/100),1);

pic_num=1;
set(handles.Stop_button,'userdata',0);
for k=1:nt
    if get(handles.Stop_button, 'userdata') % stop condition
		break;
	end
    t=k*dt;
    c_p=[c(2:end) c(1)];  % c(i+1)
    c_l=[c(end) c(1:end-1)]; % c(i-1)
    
    if eqn_type == 1 % A-C equation
        c=c-dt*M*[df*(4*c.^3-4*c)-gamma*(c_p-2*c+c_l)/dx/dx]; % update c by the Allen-Cahn equation
    elseif eqn_type == 2 % C-H equation
        miu = df*(4*c.^3-4*c)-gamma*(c_p-2*c+c_l)/dx/dx;
        miu_p = [miu(2:end) miu(1)];  % miu(i+1)
        miu_l = [miu(end) miu(1:end-1)]; % miu(i-1)
        c=c+dt*M*(miu_p-2*miu+miu_l)/dx/dx;
    end
    c(c>0.9999)=0.9999; c(c<-0.9999)=-0.9999;
    if mod(k,nprint)==0 % plot the results and save animation
        if get(handles.Stop_button, 'userdata') % stop condition
            break;
        end
        plot(x,c0,'k--',x,c,'k-','linewidth',1.0);
        legend('Initial state','Instant state');
        title(sprintf('Time =%1.2e sec', t),'position',[0.5*L,max(C1,C2)+0.3*abs(C1-C2)]);
        xlabel('x');
        ylabel('\phi');
        set(gca,'TickDir','out');
        axis([0 L min(C1,C2)-abs(C1-C2)*0.5 max(C1,C2)+abs(C1-C2)*0.5]);
        set(gca,'Fontname','Times New Roman','Fontsize',17);
        drawnow;
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        if pic_num == 1
            imwrite(I,map,'result.gif','gif', 'Loopcount',inf,'DelayTime',0.);
        else
            imwrite(I,map,'result.gif','gif','WriteMode','append','DelayTime',0.);
        end
        pic_num = pic_num + 1;
    end
end


% --- Executes on button press in Stop_button.
function Stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to Stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Stop_button,'userdata',1);
