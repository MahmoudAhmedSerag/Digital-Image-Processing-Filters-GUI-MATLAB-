function varargout = Brain_Last(varargin)
% BRAIN_LAST MATLAB code for Brain_Last.fig
%      BRAIN_LAST, by itself, creates a new BRAIN_LAST average raises the existing
%      singleton*.
%
%      H = BRAIN_LAST returns the handle to a new BRAIN_LAST average the handle to
%      the existing singleton*.
%
%      BRAIN_LAST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAIN_LAST.M with the given input arguments.
%
%      BRAIN_LAST('Property','Value',...) creates a new BRAIN_LAST average raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Brain_Last_OpeningFcn gets called.  An
%      unrecognized property name average invalid value makes property application
%      stop.  All inputs are passed to Brain_Last_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Brain_Last

% Last Modified by GUIDE v2.5 19-May-2022 12:11:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Brain_Last_OpeningFcn, ...
                   'gui_OutputFcn',  @Brain_Last_OutputFcn, ...
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


% --- Executes just before Brain_Last is made visible.
function Brain_Last_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Brain_Last (see VARARGIN)

% Choose default command line output for Brain_Last
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Brain_Last wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Brain_Last_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Median.
function Median_Callback(~, ~, handles)
% hObject    handle to Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1 

[path, nofile] = imgetfile();

if nofile
  msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end

img1 = imread(path);
axes(handles.axes1);
imshow(img1);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');

axes (handles.axes2)
if size(img1,3) == 3
    img1 = rgb2gray(img1);
end

K = medfilt2(img1);
axes(handles.axes2);
imshow(K);
title('\fontsize{10}\color[rgb]{1,1,1} Med Filtter');
imsave



% --- Executes on button press in Negative.
function Negative_Callback(~, ~, handles)
% hObject    handle to Negative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global img3 img4

[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end

img3 = imread(path);
img4 = img3;

axes(handles.axes9);
imshow(img3);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image'); 
% levels of the 8-bit image
L = 2 ^ 8;     
% finding the negative                   
img3 = (L - 1) - img4;
axes(handles.axes3);
imshow(img3);
title('\fontsize{10}\color[rgb]{1,1,1} Negative Image');
imsave



% --- Executes on button press in Histogram.
function Histogram_Callback(~, ~, handles)
% hObject    handle to Histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img1 

[path, nofile] = imgetfile();

if nofile
  msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end

img1 = imread(path);

axes(handles.axes4);
imshow(img1);
axes(handles.axes5);
imhist(img1)

    
% --- Executes on button press in Gray_Level.
function Gray_Level_Callback(~, ~, handles)
% hObject    handle to Gray_Level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global img9 

[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end

img9 = imread(path);
axes(handles.axes8);
imshow(img9);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');

gray_image = double(img9);
[rows,cols]=size(gray_image);
mask = [0,1,0;1,-4,1;0,1,0];
img9 = gray_image;
for i=2:rows-1
 for j=2:cols-1
     temp = mask.*gray_image(i-1:i+1,j-1:j+1);
     value = sum(temp(:));
     img9(i, j)= value;
end
end
img9 = uint8(img9);
axes(handles.axes12)
imshow(img9);
title('\fontsize{10}\color[rgb]{1,1,1} sharp image');
imsave



function Edge_Detection_Callback(~, ~, handles) %#ok<*DEFNU>
% hObject    handle to Edge_Detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img7 img8
[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end
img7 = imread(path);
img7 = im2double(img7);
img8 = img7;

axes(handles.axes6);
imshow(img7);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');

if size(img7,3) == 3
    img7 = rgb2gray(img7);
end
K = medfilt2(img7);
C = double(K);

for i = 1:size(C,1)-2
    for j = 1:size(C,2)-2
        %Sobel mask for X-direction
        Gx = ((2*C(i+2,j+1)+ C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Sobel mask for Y-Dirction
        Gy = ((2*C(i+1,j+2)+ C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
        
        %The Gradient of The Image
        %B(i,j) abs(Gx) + abs(Gy)
        B(i,j) = sqrt(Gx.^2+Gy.^2);
    end
end
axes(handles.axes11)
imshow(B);
title('\fontsize{10}\color[rgb]{1,1,1} Edge Detection');
imsave


% --- Executes on button press in stretch.
function stretch_Callback(~, ~, handles)
% hObject    handle to stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img11 
[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end
img11 = imread(path);

axes(handles.axes13);
imshow(img11);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');
stretched_Image = imadjust(img11, stretchlim(img11, [0.05, 0.95]),[]);

%subplot(2,2,1), imshow(image), title('Original Image');
%subplot(2,2,2), 
axes(handles.axes14);
imshow(stretched_Image),
title('\fontsize{10}\color[rgb]{1,1,1} Strethced Image');
imsave


% --- Executes on button press in Log.
function Log_Callback(~, ~, handles)
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img13 
[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end
img13 = imread(path);

axes(handles.axes15);
imshow(img13);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');

b=im2double(img13);
s=(2*log(1+b))*256;
s1=uint8(s);
axes(handles.axes16);
imshow(s1);
title ('\fontsize{10}\color[rgb]{1,1,1} c=2');
imsave

%sp=(2*log(1+b))*256;
%s2=uint8(sp)
%subplot(2,2,3)
%imshow(s2);
%title 'c=2'


%sp2=(3*log(1+b))*256;
%s3=uint8(sp2)
%subplot(2,2,4)
%imshow(s3);
%title 'c=3'


% --- Executes on button press in thresholding.
function thresholding_Callback(~, ~, handles)
% hObject    handle to thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%x = imread('grayImage.png');
global img15 
[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end
img15 = imread(path);

axes(handles.axes17);
imshow(img15);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');
y=img15;
[w, h]=size(img15);
for i=1:w
    for j=1:h
        if img15(i,j)>=150
            y(i,j)=255;
        else
            y(i,j)=0;
        end    
    end
end

axes(handles.axes18);
imshow(y);
title('\fontsize{10}\color[rgb]{1,1,1} Thresholding');
imsave


% --- Executes on button press in Down_Sample.
function Down_Sample_Callback(~, ~, handles)
% hObject    handle to Down_Sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img17 
[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end
img17 = imread(path);
%img18=rgb2gray(img17);
if size(img17,3) == 3
    img17 = rgb2gray(img17);
end
axes(handles.axes19);
imshow(img17);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');
[r,c]=size(img17);
%Aa=A(1:2:r,1:2:c)
%figure,imshow(Aa);
%title('downSampled Image 512*512');
%Aaa=A(1:4:r,1:4:c);
%figure,imshow(Aaa);
%title('downsampled Image 256*256');
Aaa=img17(1:8:r,1:8:c);
axes(handles.axes20);
imshow(Aaa);
title('\fontsize{10}\color[rgb]{1,1,1} downsampled Image 128*128');
imsave


% --- Executes on button press in Min_Max.
function Min_Max_Callback(~, ~, handles)
% hObject    handle to Min_Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%original =imread('mf.png')
global img19 
[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end
img19 = imread(path);

BW = im2bw(img19,0.6);           %Convert into black and white image

minf=@(x) min(x(:));                %set 'min()' filter
maxf=@(x)max(x(:));                 %set 'max()' filter
min_Image=nlfilter(BW,[3 3],minf);  %Apply over 3 x 3 neighbourhood
max_Image=nlfilter(BW,[3 3],maxf);  %Apply over 3 x 3 neighbourhood
axes(handles.axes21);
imshow(BW), 
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');   %Display image 
axes(handles.axes22);
imshow(min_Image), 
title('\fontsize{10}\color[rgb]{1,1,1} Min'); %Display min image
axes(handles.axes23);
imshow(max_Image),
title('\fontsize{10}\color[rgb]{1,1,1} Max'); %Display max image
imsave


function Gray_Level_Slicing_Callback(~, ~, handles)
% hObject    handle to Gray_Level_Slicing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img5 img6

[path, nofile] = imgetfile();

if nofile
    msgbox(sprintf('Image not found!!!'), 'Error', 'Warning');
    return
end

img5 = imread(path);
img5 = im2double(img5);
img6 = img5;

axes(handles.axes7);
imshow(img5);
title('\fontsize{10}\color[rgb]{1,1,1} Original Image');

img5 = imread(path);
%img5 = rgb2gray(img5);  
if size(img5,3) == 3
    img5 = rgb2gray(img5);
end
x = img5;
[rows, cols] = size(img5);
for row_index=1:1:rows
    for col_index=1:1:cols
        if(img5(row_index,col_index)>=100 && img5(row_index,col_index)<=150)
            x(row_index,col_index) = 255;
        else
             x(row_index,col_index) = 0;
        end
    end
end
axes(handles.axes10)
imshow(x);
title('\fontsize{10}\color[rgb]{1,1,1} Gray Level Slicing Image');
imsave
