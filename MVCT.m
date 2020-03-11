
function varargout = MVCT(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MVCT_OpeningFcn, ...
    'gui_OutputFcn',  @MVCT_OutputFcn, ...
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

end

% --- Executes just before MVCT is made visible.
function MVCT_OpeningFcn(hObject, ~, handles, varargin)
global myDocsFolder;
global refdataAB;
global refdataBC;
global refdataCA;
global refdataBD;
global refdataMean;
global refdataStdev;
global refdataMaterial;
global currentFolder;
userProfile = getenv('USERPROFILE');
myDocsFolder = sprintf('%s\\Documents\\TomoMQA\\', userProfile);
mkdir(myDocsFolder);
currentFolder=0;
refdataAB=0;
refdataBC=0;
refdataCA=0;
refdataBD=0;
refdataMean=zeros(6,1);
refdataStdev=zeros(6,1);
refdataMaterial=zeros(20,4);

% check qaparam.dat is exist
if isfile([myDocsFolder 'qaparam.dat'])
    % nothing do
else
    fid=fopen([myDocsFolder 'qaparam.dat'], 'w');
    fprintf(fid, '5\n120\n70\n0\n1\n30\n50\n1400\n');    
    fclose(fid);
end

% Choose default command line output for MVCT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% Update axes, table, and radiobox for initialization
windows_initial(handles);
visible_off(handles);

end %MVCT_OpeningFcn

% --- Outputs from this function are returned to the command line.
function varargout = MVCT_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end %varargout

function [threedarray, xthickness, ythickness, zthickness, rescaleIntercept] = gatherImages(folder)
global acquisitionDate;
%GATHERIMAGES looks through a folder, gets all DICOM files and assembles
%them into a viewable 3d format.

d = sortDirectory(folder); %Sort in ascending order of instance number
topimage = dicomread(d(1,:));
metadata = dicominfo(d(1,:));
acquisitionDate= metadata.AcquisitionDate;
[group1, element1] = dicomlookup('PixelSpacing');
[group2, element2] = dicomlookup('SliceThickness');
resolution = metadata.(dicomlookup(group1, element1));
xthickness = resolution(1); ythickness = resolution(2);
zthickness = metadata.(dicomlookup(group2, element2));
rescaleIntercept=metadata.RescaleIntercept;

threedarray = zeros(size(topimage,1),size(topimage,2),size(d,1));
threedarray(:,:,1) = topimage;

for i = 2:size(d,1)
    threedarray(:,:,i) = dicomread(d(i,:));
end

end %end gatherImages(folder)

function d = sortDirectory(folder)
%SORTDIRECTORY sorts based on instance number

cd(folder);
d = ls('*.dcm');
m = size(d,1);

[group, element] = dicomlookup('InstanceNumber');
sdata(m) = struct('imagename','','instance',0);

for i = 1:m
    metadata = dicominfo(d(i,:));
    position = metadata.(dicomlookup(group, element));
    sdata(i) = struct('imagename',d(i,:),'instance',position);
end

[~, order] = sort([sdata(:).instance],'ascend');
sorted = sdata(order).';

for i = 1:m
    d(i,:) = sorted(i).imagename;
end

end %end sortDirectory(dir)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% call global variables
global myDocsFolder;
global threedarray;
global xthickness;
global ythickness;
global zthickness;
global mPoint;
global startR;
global maxHU;
global fiducialR;
global rescaleIntercept;
global materialData;
global currentFolder;
global step2_flag;
step2_flag=0;
global refdataAB;
global refdataBC;
global refdataCA;
global refdataBD;
global refdataMean;
global refdataStdev;
global refdataMaterial;
global refacquisitionDate;
global spatialSliceNo;

set(handles.text16,'String','Reference data is not loaded.');

refdataAB=0;
refdataBC=0;
refdataCA=0;
refdataBD=0;
refdataMean=zeros(6,1);
refdataStdev=zeros(6,1);
refdataMaterial=zeros(20,4);
refacquisitionDate='00000000';

windows_initial(handles);
visible_off(handles);
% change mouse pointer
set(handles.figure1, 'pointer', 'arrow');
oldpointer = get(handles.figure1, 'pointer');
set(handles.figure1, 'pointer', 'watch');

try
    if(currentFolder==0) 
        folder = uigetdir(myDocsFolder);
    else
        folder = uigetdir(currentFolder);
    end
    currentFolder=folder;
    w = cd; cd(folder);
catch %#ok<CTCH>
    return
end
try
    [threedarray, xthickness, ythickness, zthickness, rescaleIntercept] = gatherImages(folder);
catch %#ok<CTCH>
    msgbox('You selected a folder that does not contain .dcm images');
    set(handles.figure1, 'pointer', 'arrow');
    threedarray = NaN; xthickness = NaN; ythickness = NaN; zthickness = NaN;
    cd(w);
    return
end
cd(w);
[~, Q, R] = size(threedarray);

% Boundary Check of Cheese Phantom
checkPoint1 = [256-40 256-40];
checkPoint2 = [256+40 256-40];
checkPoint3 = [256-40 256+40];
checkPoint4 = [256+40 256+40];
flag=0;
startR = 1;

for i=1:R
    checkHU1= threedarray(checkPoint1(1), checkPoint1(2), i);
    checkHU2= threedarray(checkPoint2(1), checkPoint2(2), i);
    checkHU3= threedarray(checkPoint3(1), checkPoint3(2), i);
    checkHU4= threedarray(checkPoint4(1), checkPoint4(2), i);
    if(flag==0) % start point
        if(checkHU1 > 200  && checkHU2 > 200 && checkHU3 > 200 && checkHU4 > 200 )
            flag=1;
            startR = i;
        end
    end
end

% Find fiducial markers
[maxHU, fiducialR] = max(max(max(threedarray(1:150, :, :))));
[maxHU2, fiducialR2] = max(max(max(threedarray(end-150:end, :, :))));
[maxHU3, fiducialR3] = max(max(max(threedarray(:, 1:150, :))));
[maxHU4, fiducialR4] = max(max(max(threedarray(:, end-150:end, :))));


if(maxHU<maxHU2)
    maxHU = maxHU2;
    fiducialR=fiducialR2;
end

if(maxHU<maxHU3)
    maxHU = maxHU3;
    fiducialR=fiducialR3;
    
end

if(maxHU<maxHU4)
    maxHU = maxHU4;
    fiducialR=fiducialR4;
end

% --------------
% Image Rotation
% --------------
flag1 = sum(sum(sum(threedarray(1:150,:,fiducialR)>maxHU-800)));
flag2 = sum(sum(sum(threedarray(end-150:end,:,fiducialR)>maxHU-800)));
flag3 = sum(sum(sum(threedarray(:,1:150,fiducialR)>maxHU-800)));
flag4 = sum(sum(sum(threedarray(:,end-150:end,fiducialR)>maxHU-800)));
maxflag = flag1;
flag = 1;
if(maxflag<flag2)
    maxflag=flag2;
    flag = 2;
end
if(maxflag<flag3)
    maxflag=flag3;
    flag = 3;
end
if(maxflag<flag4)
    flag=4;
end
if(flag==2)
    for i=1:R
        threedarray(:,:,i) = threedarray(:,end:-1:1,i)';
        threedarray(:,:,i) = threedarray(:,end:-1:1,i)';
    end
elseif(flag==3)
    for i=1:R
        threedarray(:,:,i) = threedarray(:,:,i)';
    end
elseif(flag==4)
    for i=1:R
        threedarray(:,:,i) = threedarray(:,end:-1:1,i)';
    end
end

% -----------------------------------------------------------
% length Calculation - TG-148 VI.B.1.a. Geometric distortions
% -----------------------------------------------------------
% draw fiducial (A,B,C) - axial
maxHU_Tol = (str2double(get(handles.edit7,'String')));
while(1)
    measurements = regionprops(threedarray(:,:,fiducialR)>maxHU_Tol, 'Centroid');
    if(length(measurements)>=3 && length(measurements)<=4) 
        break;
    else
       maxHU_Tol= maxHU_Tol+50;
    end
end
mPoint=0;
idx = 1;
for k = 1:length(measurements)
    if((length(measurements)==4 && k~=2) || length(measurements)==3)
        temp = measurements(k).Centroid;
        mPoint(idx,1) = temp(1);
        mPoint(idx,2) = temp(2);
        idx=idx+1;
    end
end
mPoint(idx,:)=mPoint(1,:);

% --------------
%  axial image - axis1
% --------------
% set colormap (like CT)
colormap('Gray'); %Dealing only with grayscaled images (*.DCM)

% draw CT image
imagesc(threedarray(:,:,fiducialR),'Parent',handles.axes1);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes1);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes1);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes1);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes1);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes1);

% --------------
% coronal image - axis2
% --------------
% find coronal CT image

fiducialArray2=squeeze(threedarray(:,int16(mPoint(2,1)),:))';
slideNo=sprintf("IM: %d/%d",int16(mPoint(2,1)), Q);

% set array size
coronalImage = zeros(512,512);
[imgX, imgY] = size(fiducialArray2);
increasNo = int16(512*0.5)-int16(imgX*0.5);
coronalImage(increasNo:increasNo+imgX-1, 1:imgY) = fiducialArray2;

% draw CT image (coronal)
imagesc(coronalImage,'Parent', handles.axes2);

% direction label
text(250,15,'S','FontSize',12, 'Color','Y', 'Parent', handles.axes2);
text(250,490,'I','FontSize',12, 'Color','Y', 'Parent', handles.axes2);
text(10,250,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes2);
text(490,250,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes2);

% image slide NO
text(10,15, slideNo,'FontSize',13, 'Color','r', 'Parent', handles.axes2);


% set colormap limits (=axes1)
caxis(handles.axes2, caxis(handles.axes1));


% --------------
% Uniformity - axis3
% --------------

% draw CT image
imagesc(threedarray(:,:,fiducialR-10),'Parent',handles.axes3);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes3);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes3);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes3);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes3);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-10, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes3);

% set colormap limits (=axes1)
caxis(handles.axes3, caxis(handles.axes1));


% --------------
% Noise - axis4
% --------------

% draw CT image
imagesc(threedarray(:,:,fiducialR-10),'Parent',handles.axes4);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes4);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes4);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes4);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes4);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-10, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes4);

% set colormap limits (=axes1)
caxis(handles.axes4, caxis(handles.axes1));


% --------------
% Contrast - axes5
% --------------
spatial_image=threedarray(:,:,fiducialR-35);
if(fiducialR-95>0) 
    spatial_image2=threedarray(:,:,fiducialR-95);
end
contrast_image=threedarray(:,:,fiducialR-50);

% draw CT image
imagesc(contrast_image,'Parent',handles.axes5);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes5);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes5);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes5);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes5);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-50, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes5);

caxis(handles.axes5, caxis(handles.axes1));

% Numbering.... - axes 6
imagesc(spatial_image,'Parent',handles.axes6);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes6);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-35, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes6);

caxis(handles.axes6, caxis(handles.axes1));

radius_inner = 65;
radius_outer = 110;
plugPoint=zeros(20,2);

uniform_image=threedarray(:,:,fiducialR-10);
boundary = uniform_image>800;
hold(handles.axes3, 'on');
centers = regionprops(boundary, 'Centroid');

for k = 1:length(centers)
    if(sum(centers(k).Centroid > [200, 200]) + sum(centers(k).Centroid < [300, 300]) == 4)
        centerPoint = centers(k).Centroid;
        break;
    end
end

% Find Plug Points
for i=1:8
    plugPoint(i,:) = [ centerPoint(1,1) + cos(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/xthickness, centerPoint(1,2) + sin(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/ythickness];
end

for i=1:12
    plugPoint(8+i,:) = [centerPoint(1,1) + cos(2*pi()*i/12 - 15*pi()/180)*radius_outer/xthickness, centerPoint(1,2) + sin(2*pi()*i/12 - 15*pi()/180)*radius_outer/ythickness];
end

% Write number of plug points
for i=1:20
    if (i<10)
        minuspoint = 5;
    else
        minuspoint = 10;
    end
    text(plugPoint(i,1)-minuspoint,plugPoint(i,2),num2str(i),'FontSize',10, 'Color','W', 'Parent', handles.axes6);
end


% ------------------
% Contrast - Table1
% ------------------
materialData=zeros(20,4);

tabledataHTML = reshape(strtrim(cellstr(num2str(materialData(:),'%.2f'))), size(materialData));
tabledataHTML = strcat('<html><div style="width: 90; text-align: center;">', tabledataHTML,'</div></html>');
set(handles.uitable1, 'Data', tabledataHTML);

%caxis(handles.axes7, caxis(handles.axes1));


% ----------------------------
% Spatial Resolution - axes7
% ----------------------------


% label for image slide number
[columnsInImage, rowsInImage] = meshgrid(1:512, 1:512);
radius=uint16(str2double(get(handles.edit2,'String')));
radius_inner = 65;
radius_outer = 110;

consistancy_average = zeros(20,1);
consistancy_stdev = zeros(20,1);
consistancy_max = zeros(20,1);
consistancy_min = zeros(20,1);

% spatial resolution - inside
plugPoint=zeros(20,2);
for i=1:8
    plugPoint(i,:) = [ centerPoint(1,1) + cos(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/xthickness, centerPoint(1,2) + sin(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/ythickness];
    circlePixels = (rowsInImage - plugPoint(i,2)).^2 + (columnsInImage - plugPoint(i,1)).^2 <= (radius/xthickness).^2;
    t=spatial_image.*circlePixels;
    consistancy_average(i) = mean(t(t>0))+rescaleIntercept;
    consistancy_stdev(i) = std(t(t>0));
    consistancy_max(i) = max(t(t>0));
    consistancy_min(i) = min(t(t>0));
end

for i=1:12
    plugPoint(8+i,:) = [centerPoint(1,1) + cos(2*pi()*i/12 - 15*pi()/180)*radius_outer/xthickness, centerPoint(1,2) + sin(2*pi()*i/12 - 15*pi()/180)*radius_outer/ythickness];
    circlePixels = (rowsInImage - plugPoint(8+i,2)).^2 + (columnsInImage - plugPoint(8+i,1)).^2 <= (radius/xthickness).^2;
    t=spatial_image.*circlePixels;
    consistancy_average(8+i) = mean(t(t>0))+rescaleIntercept;
    consistancy_stdev(8+i) = std(t(t>0));
    consistancy_max(8+i) = max(t(t>0));
    consistancy_min(8+i) = min(t(t>0));
end
materialData(:,1)= consistancy_average;
materialData(:,2)= consistancy_stdev;
[~, resol_plug]=max(materialData(:,2));
resolution_plug_no  = resol_plug;
spatialSliceNo = 35;
consistancy_average = zeros(20,1);
consistancy_stdev = zeros(20,1);
consistancy_max = zeros(20,1);
consistancy_min = zeros(20,1);

% spatial resolution - outside
if(fiducialR-95>0) 
    plugPoint2=zeros(20,2);
    for i=1:8
        plugPoint2(i,:) = [ centerPoint(1,1) + cos(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/xthickness, centerPoint(1,2) + sin(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/ythickness];
        circlePixels = (rowsInImage - plugPoint2(i,2)).^2 + (columnsInImage - plugPoint2(i,1)).^2 <= (radius/xthickness).^2;
        t=spatial_image2.*circlePixels;
        consistancy_average(i) = mean(t(t>0))+rescaleIntercept;
        consistancy_stdev(i) = std(t(t>0));
        consistancy_max(i) = max(t(t>0));
        consistancy_min(i) = min(t(t>0));
    end

    for i=1:12
        plugPoint2(8+i,:) = [centerPoint(1,1) + cos(2*pi()*i/12 - 15*pi()/180)*radius_outer/xthickness, centerPoint(1,2) + sin(2*pi()*i/12 - 15*pi()/180)*radius_outer/ythickness];
        circlePixels = (rowsInImage - plugPoint2(8+i,2)).^2 + (columnsInImage - plugPoint2(8+i,1)).^2 <= (radius/xthickness).^2;
        t=spatial_image2.*circlePixels;
        consistancy_average(8+i) = mean(t(t>0))+rescaleIntercept;
        consistancy_stdev(8+i) = std(t(t>0));
        consistancy_max(8+i) = max(t(t>0));
        consistancy_min(8+i) = min(t(t>0));
    end
end
    if(max(materialData(:,2)) < max(consistancy_stdev))
        materialData(:,1)= consistancy_average;
        materialData(:,2)= consistancy_stdev;    
        [~, resol_plug]=max(materialData(:,2));
        resolution_plug_no  = resol_plug;
        spatialSliceNo = 95;

        spatial_image=spatial_image2;

        % re-drawing spatial image to Axis6

        % Numbering.... - axes 6
        imagesc(spatial_image,'Parent',handles.axes6);

        % labels for direction
        text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
        text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
        text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
        text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes6);

        % label for image slide number
        txt=sprintf("IM: %d/%d",fiducialR-spatialSliceNo, R);
        text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes6);

        caxis(handles.axes6, caxis(handles.axes1));
    end


    if(resolution_plug_no <1 && resolution_plug_no >20)
        resolution_plug_no =1;
    end
    set(handles.edit1,'String',num2str(resolution_plug_no));
    resolution_image = spatial_image(uint16(plugPoint(resolution_plug_no,2))-20:uint16(plugPoint(resolution_plug_no,2))+20,...
        uint16(plugPoint(resolution_plug_no,1))-20:uint16(plugPoint(resolution_plug_no,1))+20);
    imagesc(resolution_image,'Parent',handles.axes7);


% label for image slide number

txt=sprintf("IM: %d/%d",fiducialR-spatialSliceNo, R);
text(0.8,1.5,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes7);

%caxis(handles.axes7, caxis(handles.axes1));

% delete axis ticks
visible_off(handles);

% mouse pointer on
set(handles.figure1, 'pointer', oldpointer)

end % pushbutton1_Callback


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% call global variables
global threedarray;
global xthickness;
global ythickness;
global zthickness;
global mPoint;
global startR;
global maxHU;
global fiducialR;
global rescaleIntercept;
global materialData;
global step2_flag;

% for data
global dataAB;
global dataBC;
global dataCA;
global dataBD;
global dataMean;
global dataStdev;
global maxDiff_Uniform;
% for refdata
global refdataAB;
global refdataBC;
global refdataCA;
global refdataBD;
global refdataMean;
global refdataStdev;
global refdataMaterial;
global maxDiff_Uniform_ref;



step2_flag=1;

[~, Q, R] = size(threedarray);

% -----------------------------------------------------------
% length Calculation - TG-148 VI.B.1.a. Geometric distortions
% -----------------------------------------------------------

% draw fiducial (A,B,C,D) - axial
lengthBC = sqrt (power((mPoint(1,1)- mPoint(2,1))*xthickness,2) + power((mPoint(1,2)- mPoint(2,2))*ythickness,2));
lengthAB = sqrt (power((mPoint(2,1)- mPoint(3,1))*xthickness,2) + power((mPoint(2,2)- mPoint(3,2))*ythickness,2));
lengthCA = sqrt (power((mPoint(3,1)- mPoint(4,1))*xthickness,2) + power((mPoint(3,2)- mPoint(4,2))*ythickness,2));
lengthBD = abs(startR-fiducialR)*zthickness;


% --------------
%  axial image
% --------------
% set colormap (like CT)
colormap('Gray'); %Dealing only with grayscaled images (*.DCM)

% draw CT image
imagesc(threedarray(:,:,fiducialR),'Parent',handles.axes1);

% find fiducial boundary
C = bwboundaries(threedarray(:,:,fiducialR)>maxHU-800,'noholes');
hold (handles.axes1, 'on');

% draw boundary contours
for k = 1:length(C)
    boundary = C{k};
    plot(boundary(:,2), boundary(:,1),'r','LineWidth',1 ,'Parent', handles.axes1);
end

% draw center points of boundaries
for k = 1:3
    scatter(mPoint(k,1), mPoint(k,2), 'b','LineWidth',0.5,'marker','.', 'MarkerEdgeColor','Blue','MarkerFaceColor','Blue','Parent',handles.axes1);
end
% draw measurement lines (A-B, B-C, and C-A)
plot(mPoint(:, 1), mPoint(:, 2), 'g','LineWidth',1, 'Parent', handles.axes1);

text(mPoint(1,1)-20,mPoint(1,2),'C','FontSize',13, 'Color','w', 'Parent', handles.axes1);
text(mPoint(2,1)-10,mPoint(2,2)-20,'B','FontSize',13, 'Color','w', 'Parent', handles.axes1);
text(mPoint(3,1)+10,mPoint(3,2),'A','FontSize',13, 'Color','w', 'Parent', handles.axes1);

% label for results (A-B, B-C, and C-A)
dataAB=lengthAB/10;
txt=sprintf("A-B: %.2f cm",dataAB);
text(330,400,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes1);
dataBC=lengthBC/10;
txt=sprintf("B-C: %.2f cm",dataBC);
text(330,400+40,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes1);
dataCA=lengthCA/10;
txt=sprintf("C-A: %.2f cm",dataCA);
text(330,400+80,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes1);
% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes1);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes1);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes1);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes1);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes1);

hold(handles.axes1, 'off');

% --------------
% coronal image
% --------------

% find coronal CT image

fiducialArray2=squeeze(threedarray(:,int16(mPoint(2,1)),:))';
slideNo=sprintf("IM: %d/%d",int16(mPoint(2,1)), Q);

% set array size
coronalImage = zeros(512,512);
[imgX, imgY] = size(fiducialArray2);
increasNo = int16(512*0.5)-int16(imgX*0.5);
coronalImage(increasNo:increasNo+imgX-1, 1:imgY) = fiducialArray2;

% draw CT image (coronal)
imagesc(coronalImage,'Parent', handles.axes2);

% change level width from CT image
fiducialArray2=fiducialArray2>2000;

% find fiducial marker (B)
C = bwboundaries(fiducialArray2,'noholes');

hold(handles.axes2, 'on');

% draw boundary contour (B)
for k = 1:length(C)
    boundary = C{k};
    plot(boundary(:,2), double(increasNo)+boundary(:,1),'r','LineWidth',1 ,'Parent', handles.axes2);
end

% line drawing from B to E
plot([mPoint(2,2) mPoint(2,2)],[double(increasNo)+fiducialR double(increasNo)+startR], 'g','LineWidth',1 ,'Parent', handles.axes2);


% marker label
text(mPoint(2,2)-20,double(increasNo)+fiducialR-5,'B','FontSize',13, 'Color','w','Parent', handles.axes2);
text(mPoint(2,2)-20,double(increasNo)+startR,'D','FontSize',13, 'Color','w','Parent', handles.axes2);

% direction label
text(250,15,'S','FontSize',12, 'Color','Y', 'Parent', handles.axes2);
text(250,490,'I','FontSize',12, 'Color','Y', 'Parent', handles.axes2);
text(10,250,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes2);
text(490,250,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes2);

% image slide NO
text(10,15,slideNo,'FontSize',13, 'Color','r', 'Parent', handles.axes2);

% label for result (B-D)
dataBD=lengthBD/10;
txt=sprintf("B-D: %.2f cm",dataBD);
text(340,360+120,txt,'FontSize',13, 'Color','r','Parent', handles.axes2);

hold (handles.axes2, 'off');

% -------------------
% Check Consistancy - Geometric Distortions
% -------------------
% pass = 1, fail = 0
gd_passfail = 1;

geometric_limit=0.2;
if(get(handles.radiobutton11, 'Value'))
    geometric_limit=0.1;
end

gd_passfail= gd_passfail && (abs(dataAB-refdataAB)-geometric_limit)<1e-10;
gd_passfail= gd_passfail && (abs(dataBC-refdataBC)-geometric_limit)<1e-10;
gd_passfail= gd_passfail && (abs(dataCA-refdataCA)-geometric_limit)<1e-10;
gd_passfail= gd_passfail && (abs(dataBD-refdataBD)-geometric_limit)<1e-10;
if(gd_passfail)     
    set(handles.radiobutton1,'Value',1);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton2,'FontWeight','normal');
    set(handles.radiobutton1,'FontWeight','bold');
else
    set(handles.radiobutton1,'Value',0);
    set(handles.radiobutton2,'Value',1);
    set(handles.radiobutton1,'FontWeight','normal');
    set(handles.radiobutton2,'FontWeight','bold');
end

% --------------
% Uniformity - Axis 3, 4
% --------------

% draw CT image
uniform_image=threedarray(:,:,fiducialR-10);
imagesc(uniform_image,'Parent',handles.axes3);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes3);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes3);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes3);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes3);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-10, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes3);


boundary = uniform_image>800;
hold(handles.axes3, 'on');
centers = regionprops(boundary, 'Centroid');

for k = 1:length(centers)
    if(sum(centers(k).Centroid > [200, 200]) + sum(centers(k).Centroid < [300, 300]) == 4)
        centerPoint = centers(k).Centroid;
        break;
    end
end

radius_big=uint16(str2double(get(handles.edit3,'String')));

centerPoint(2,:) = centerPoint(1)+[0 radius_big/xthickness];
centerPoint(3,:) = centerPoint(1)-[0 radius_big/xthickness];
centerPoint(4,:) = centerPoint(1)+[radius_big/xthickness 0];
centerPoint(5,:) = centerPoint(1)-[radius_big/xthickness 0];

radius=uint16(str2double(get(handles.edit2,'String')));

viscircles(handles.axes3, centerPoint(1,:), radius/xthickness, 'linewidth',0.3);
viscircles(handles.axes3, centerPoint(2,:), radius/xthickness, 'linewidth',0.3);
viscircles(handles.axes3, centerPoint(3,:), radius/xthickness, 'linewidth',0.3);
viscircles(handles.axes3, centerPoint(4,:), radius/xthickness, 'linewidth',0.3);
viscircles(handles.axes3, centerPoint(5,:), radius/xthickness, 'linewidth',0.3);


[columnsInImage, rowsInImage] = meshgrid(1:512, 1:512);
circlePixels = (rowsInImage - centerPoint(1,2)).^2 + (columnsInImage - centerPoint(1,1)).^2 <= (radius/xthickness).^2;
t=uniform_image.*circlePixels;
uniform_average(1) = mean(t(t>0))+rescaleIntercept;
uniform_stdev(1) = std(t(t>0));

circlePixels = (rowsInImage - centerPoint(2,2)).^2 + (columnsInImage - centerPoint(3,1)).^2 <= (radius/xthickness).^2;
t=uniform_image.*circlePixels;
uniform_average(2) = mean(t(t>0))+rescaleIntercept;
uniform_stdev(2) = std(t(t>0));

circlePixels = (rowsInImage - centerPoint(3,2)).^2 + (columnsInImage - centerPoint(3,1)).^2 <=  (radius/xthickness).^2;
t=uniform_image.*circlePixels;
uniform_average(3) = mean(t(t>0))+rescaleIntercept;
uniform_stdev(3) = std(t(t>0));

circlePixels = (rowsInImage - centerPoint(4,2)).^2 + (columnsInImage - centerPoint(4,1)).^2 <= (radius/xthickness).^2;
t=uniform_image.*circlePixels;
uniform_average(4) = mean(t(t>0))+rescaleIntercept;
uniform_stdev(4) = std(t(t>0));

circlePixels = (rowsInImage - centerPoint(5,2)).^2 + (columnsInImage - centerPoint(5,1)).^2 <= (radius/xthickness).^2;
t=uniform_image.*circlePixels;
uniform_average(5) = mean(t(t>0))+rescaleIntercept;
uniform_stdev(5) = std(t(t>0));

dataMean(1:5,1) = uniform_average; 
dataStdev(1:5,1) = uniform_stdev;
txt=sprintf("Mean: %.2f\nStd: %.2f",uniform_average(1),uniform_stdev(1));
text(centerPoint(1,1)-40,centerPoint(1,2)+40,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes3);
txt=sprintf("Mean: %.2f\nStd: %.2f",uniform_average(2),uniform_stdev(2));
text(centerPoint(2,1)-40,centerPoint(2,2)+40,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes3);
txt=sprintf("Mean: %.2f\nStd: %.2f",uniform_average(3),uniform_stdev(3));
text(centerPoint(3,1)-40,centerPoint(3,2)+40,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes3);
txt=sprintf("Mean: %.2f\nStd: %.2f",uniform_average(4),uniform_stdev(4));
text(centerPoint(4,1)-40,centerPoint(4,2)+40,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes3);
txt=sprintf("Mean: %.2f\nStd: %.2f",uniform_average(5),uniform_stdev(5));
text(centerPoint(5,1)-40,centerPoint(5,2)+40,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes3);


% set colormap limits (=axes1)
caxis(handles.axes3, caxis(handles.axes1));

% --------------
% Noise Test
% --------------

% draw CT image
imagesc(uniform_image,'Parent',handles.axes4);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes4);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes4);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes4);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes4);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-10, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes4);

hold(handles.axes4, 'on');

viscircles(handles.axes4, centerPoint(1,:), radius/xthickness, 'linewidth',0.3);
viscircles(handles.axes4, centerPoint(1,:), radius_big/xthickness, 'linewidth',0.3);

[columnsInImage, rowsInImage] = meshgrid(1:512,1:512);
circlePixels = (rowsInImage - centerPoint(1,2)).^2 + (columnsInImage - centerPoint(1,1)).^2 <= (radius/xthickness).^2;
t=uniform_image.*circlePixels;
noise_average(1) = mean(t(t>0))+rescaleIntercept;
noise_stdev(1) = std(t(t>0));

circlePixels = (rowsInImage - centerPoint(1,2)).^2 + (columnsInImage - centerPoint(1,1)).^2 <= (radius_big/xthickness).^2;
t=uniform_image.*circlePixels;
noise_average(2) = mean(t(t>0))+rescaleIntercept;
noise_stdev(2) = std(t(t>0));

dataMean(6,1) = noise_average(2); 
dataStdev(6,1) = noise_stdev(2);

txt=sprintf("Mean: %.2f\nStd: %.2f",noise_average(1),noise_stdev(1));
text(centerPoint(1,1)-40,centerPoint(1,2)+40,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes4);
txt=sprintf("Mean: %.2f\nStd: %.2f",noise_average(2),noise_stdev(2));
text(centerPoint(1,1)-40,centerPoint(1,2)-120,txt,'FontSize',11, 'Color','w', 'Parent', handles.axes4);

% set colormap limits (=axes1)
caxis(handles.axes4, caxis(handles.axes1));

% -------------------
% Check Consistancy - Uniformity and Noise
% -------------------
% pass = 1, fail = 0
uniform_passfail_forDC = 1;


maxDiff_Uniform_ref = abs(refdataMean(1) - refdataMean(2));
if(maxDiff_Uniform_ref< abs(refdataMean(1) - refdataMean(3)))
    maxDiff_Uniform_ref = abs(refdataMean(1) - refdataMean(3));
end
if(maxDiff_Uniform_ref< abs(refdataMean(1) - refdataMean(4)))
    maxDiff_Uniform_ref = abs(refdataMean(1) - refdataMean(4));
end
if(maxDiff_Uniform_ref< abs(dataMean(1) - refdataMean(5)))
    maxDiff_Uniform_ref = abs(dataMean(1) - refdataMean(5));
end



maxDiff_Uniform = abs(dataMean(1) - dataMean(2));
if(maxDiff_Uniform< abs(dataMean(1) - dataMean(3)))
    maxDiff_Uniform = abs(dataMean(1) - dataMean(3));
end
if(maxDiff_Uniform< abs(dataMean(1) - dataMean(4)))
    maxDiff_Uniform = abs(dataMean(1) - dataMean(4));
end
if(maxDiff_Uniform< abs(dataMean(1) - dataMean(5)))
    maxDiff_Uniform = abs(dataMean(1) - dataMean(5));
end

uniform_passfail_forDC = uniform_passfail_forDC && (maxDiff_Uniform <= 25);

HUDiff_other = (str2double(get(handles.edit6,'String')));
uniform_passfail_forNonDC =  1;
uniform_passfail_forNonDC = uniform_passfail_forNonDC && (maxDiff_Uniform <= HUDiff_other);

global maxDiff_Noise_Ref;
global maxDiff_Noise_BigCircle_Ref;

global maxDiff_Noise;
global maxDiff_Noise_BigCircle;
global spatialSliceNo;
maxDiff_Noise_Ref = (refdataStdev(1));

% for i=2:5
%     if(maxDiff_Noise_Ref < refdataStdev(i))
%         maxDiff_Noise_Ref = refdataStdev(i);
%     end
% end

maxDiff_Noise_BigCircle_Ref = (refdataStdev(6));



maxDiff_Noise = (dataStdev(1));
% for i=2:5
%     if(maxDiff_Noise < dataStdev(i))
%         maxDiff_Noise = dataStdev(i);
%     end
% end
maxDiff_Noise_BigCircle = (dataStdev(6));
HUDiff_noise = (str2double(get(handles.edit4,'String')));

noise_passfail = 1;
noise_passfail = noise_passfail && (abs(maxDiff_Noise) < HUDiff_noise);
noise_passfail = noise_passfail && (abs(maxDiff_Noise_BigCircle) < HUDiff_noise);
% noise_passfail = noise_passfail && (maxDiff_Noise <= maxDiff_Noise_Ref*(1+ConsistancyUncertainty));
% noise_passfail = noise_passfail && (maxDiff_Noise_BigCircle <= maxDiff_Noise_BigCircle_Ref*(1+ConsistancyUncertainty));

uniform_passfail=uniform_passfail_forNonDC;
if(get(handles.radiobutton9, 'Value'))
    uniform_passfail=uniform_passfail_forDC;
end

if(uniform_passfail && noise_passfail)
    set(handles.radiobutton3,'Value',1);
    set(handles.radiobutton4,'Value',0);
    set(handles.radiobutton4,'FontWeight','normal');
    set(handles.radiobutton3,'FontWeight','bold');
else
    set(handles.radiobutton3,'Value',0);
    set(handles.radiobutton4,'Value',1);
    set(handles.radiobutton3,'FontWeight','normal');
    set(handles.radiobutton4,'FontWeight','bold');
end

% --------------
% Contrast - axes5
% --------------
contrast_image=threedarray(:,:,fiducialR-50);
spatial_image=threedarray(:,:,fiducialR-spatialSliceNo);

% draw CT image
imagesc(contrast_image,'Parent',handles.axes5);

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes5);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes5);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes5);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes5);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-50, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes5);
hold(handles.axes5, 'on');

radius=uint16(str2double(get(handles.edit2,'String')));
radius_inner = 65;
radius_outer = 110;

consistancy_average = zeros(20,1);
consistancy_stdev = zeros(20,1);
consistancy_max = zeros(20,1);
consistancy_min = zeros(20,1);

% inner cycle
plugPoint=zeros(20,2);
for i=1:8
    plugPoint(i,:) = [ centerPoint(1,1) + cos(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/xthickness, centerPoint(1,2) + sin(2*pi()*i/8 - 22.5*pi()/180)*radius_inner/ythickness];
    viscircles(handles.axes5, plugPoint(i,:), radius/xthickness, 'linewidth',0.3);
    circlePixels = (rowsInImage - plugPoint(i,2)).^2 + (columnsInImage - plugPoint(i,1)).^2 <= (radius/xthickness).^2;
    t=contrast_image.*circlePixels;
    consistancy_average(i) = mean(t(t>0))+rescaleIntercept;
    consistancy_stdev(i) = std(t(t>0));
    consistancy_max(i) = max(t(t>0));
    consistancy_min(i) = min(t(t>0));
end

for i=1:12
    plugPoint(8+i,:) = [centerPoint(1,1) + cos(2*pi()*i/12 - 15*pi()/180)*radius_outer/xthickness, centerPoint(1,2) + sin(2*pi()*i/12 - 15*pi()/180)*radius_outer/ythickness];
    viscircles(handles.axes5, plugPoint(8+i,:), radius/xthickness, 'linewidth',0.3);
    circlePixels = (rowsInImage - plugPoint(8+i,2)).^2 + (columnsInImage - plugPoint(8+i,1)).^2 <= (radius/xthickness).^2;
    t=contrast_image.*circlePixels;
    consistancy_average(8+i) = mean(t(t>0))+rescaleIntercept;
    consistancy_stdev(8+i) = std(t(t>0));
    consistancy_max(8+i) = max(t(t>0));
    consistancy_min(8+i) = min(t(t>0));
end


caxis(handles.axes5, caxis(handles.axes1));


imagesc(spatial_image,'Parent',handles.axes6);
for i=1:20
    if (i<10)
        minuspoint = 5;
    else
        minuspoint = 10;
    end
    text(plugPoint(i,1)-minuspoint,plugPoint(i,2),num2str(i),'FontSize',10, 'Color','W', 'Parent', handles.axes6);
end

materialData(:,1) = consistancy_average;
materialData(:,2) = consistancy_stdev;

% -------------------
% Check Consistancy - Contrast Test
% -------------------
% pass = 1, fail = 0
contrastFlag=zeros(20,1);
contrast_passfail=1;
MaxHU_DiffForWater = uint16(str2double(get(handles.edit5,'String')));
MaxHU_DiffForOther = uint16(str2double(get(handles.edit6,'String')));
contrastFailList=zeros(20,1);
if(get(handles.radiobutton10,'value'))     
    MaxHU_DiffForWater=MaxHU_DiffForOther;
end

for i=1:20
    eval(['checkFlag = get(handles.checkbox' num2str(i) ', ''Value'');']);
    if(checkFlag==1)
        minHU_Diff=999;
        tempFlag = 0;
        for j=1:20        
            if(refdataMaterial(j)==0) 
                break;
            end
            if(minHU_Diff > abs(materialData(i)-refdataMaterial(j)) && contrastFlag(j)==0 ) 
                minHU_Diff = abs(materialData(i)-refdataMaterial(j));
                tempFlag = j;
            end
        end
        if (tempFlag ==0)  
            tempFlag = 1;
        end
        materialData(i,3) = refdataMaterial(tempFlag);
        materialData(i,4) = minHU_Diff;
        contrastFlag(tempFlag)=1;
        
        if(minHU_Diff>MaxHU_DiffForOther)
            contrast_passfail=0;
            contrastFailList(i)=1;
        end
        if(minHU_Diff>MaxHU_DiffForWater && abs(refdataMaterial(tempFlag))<30)
            contrast_passfail=0;
            contrastFailList(i)=1;
        end
        if(refdataMaterial(tempFlag) == 0)
            contrast_passfail=0;
            contrastFailList(i)=1;
        end
    else
        materialData(i,3) = 0;
        materialData(i,4) = 0;        
    end
end

if(contrast_passfail)     
    set(handles.radiobutton5,'Value',1);
    set(handles.radiobutton6,'Value',0);   
    set(handles.radiobutton6,'FontWeight','normal');
    set(handles.radiobutton5,'FontWeight','bold');
else
    set(handles.radiobutton5,'Value',0);
    set(handles.radiobutton6,'Value',1);
    set(handles.radiobutton5,'FontWeight','normal');
    set(handles.radiobutton6,'FontWeight','bold');
end

idx_fail = (contrastFailList == 1);
idx_pass = (contrastFailList == 0);
materialDataTable = reshape(strtrim(cellstr(num2str(materialData(:),'%.2f'))), size(materialData));
materialDataTable(idx_fail,:) = strcat('<html><div style="width: 90;text-align: center; color: red;">', materialDataTable(idx_fail,:),'</div></html>');
materialDataTable(idx_pass,:) = strcat('<html><div style="width: 90; text-align: center;">', materialDataTable(idx_pass,:),'</div></html>');


set(handles.uitable1, 'Data', materialDataTable);
% set colormap limits (=axes1)
caxis(handles.axes6, caxis(handles.axes1));

% labels for direction
text(250,15,'A','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
text(250,490,'P','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
text(10,250,'R','FontSize',12, 'Color','Y', 'Parent', handles.axes6);
text(490,250,'L','FontSize',12, 'Color','Y', 'Parent', handles.axes6);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-spatialSliceNo, R);
text(10,15,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes6);


resolution_plug_no  = uint8(str2double(get(handles.edit1,'String')));
if(resolution_plug_no <1 && resolution_plug_no >20)
    resolution_plug_no =1;
end

resolution_image = spatial_image(uint16(plugPoint(resolution_plug_no,2))-20:uint16(plugPoint(resolution_plug_no,2))+20,...
    uint16(plugPoint(resolution_plug_no,1))-20:uint16(plugPoint(resolution_plug_no,1))+20);
imagesc(resolution_image,'Parent',handles.axes7);

% label for image slide number
txt=sprintf("IM: %d/%d",fiducialR-spatialSliceNo, R);
text(0.8,1.5,txt,'FontSize',13, 'Color','r', 'Parent', handles.axes7);
% delete axis ticks
visible_off(handles)
end %pushbutton2_Callback

function visible_off(handles)
handles.axes1.XAxis.Visible='off';
handles.axes1.YAxis.Visible='off';
handles.axes2.XAxis.Visible='off';
handles.axes2.YAxis.Visible='off';
handles.axes3.XAxis.Visible='off';
handles.axes3.YAxis.Visible='off';
handles.axes4.XAxis.Visible='off';
handles.axes4.YAxis.Visible='off';
handles.axes5.XAxis.Visible='off';
handles.axes5.YAxis.Visible='off';
handles.axes6.XAxis.Visible='off';
handles.axes6.YAxis.Visible='off';
handles.axes7.XAxis.Visible='off';
handles.axes7.YAxis.Visible='off';
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(~, ~, handles)

% for data
global dataAB;
global dataBC;
global dataCA;
global dataBD;
global maxDiff_Uniform;
global maxDiff_Noise;
global maxDiff_Noise_BigCircle;

% for refdata
global refdataAB;
global refdataBC;
global refdataCA;
global refdataBD;
global maxDiff_Uniform_ref;
global maxDiff_Noise_Ref;
global maxDiff_Noise_BigCircle_Ref;


% Acquisition date information
global refacquisitionDate;
global acquisitionDate;
global spatialSliceNo;

set(handles.figure1, 'pointer', 'arrow');
oldpointer = get(handles.figure1, 'pointer');
set(handles.figure1, 'pointer', 'watch');

% ----- Report Generator ----- 
import mlreportgen.report.*
import mlreportgen.dom.*

rpt = Report('MVCT_Report','pdf');

tp = TitlePage;
tp.Title = 'TomoTherapy MVCT QA Report';
tp.Author = '';
tp.PubDate = datestr(now,'yyyy-mm-dd');
add(rpt,tp);

text0 = Paragraph(' ');
text0.Style = {FontSize('11')};
text0.WhiteSpace='preserve';

% Print Acquisition Date

add(rpt,text0);

bodyContent1 = {'Acquisition Date for Reference DICOM Data', [refacquisitionDate(1:4) '-' refacquisitionDate(5:6) '-' refacquisitionDate(7:8)] };
bodyContent2 = {'Acquisition Date for Current DICOM Data ' [acquisitionDate(1:4) '-' acquisitionDate(5:6) '-' acquisitionDate(7:8)]  };
emptyContent = {[],[]};
ca = [bodyContent1;bodyContent2;emptyContent];
table0 = Table(ca);
table0.Style = {FontSize('13')};

table0.Width = '100%';
table0.Children(1,1).Children(1,2).Style = {HAlign('right')};
table0.Children(1,2).Children(1,2).Style = {HAlign('right')};
add(rpt,table0);


% Print Geometric Distortion 
if(get(handles.radiobutton1,'Value'))
    gd_pass = 'Pass';
else
    gd_pass = 'Fail';
end

use_srs = get(handles.radiobutton11,'Value');

limitlength = '0.1';
if(use_srs)
    Title={'I. Geometric Distortion',  '(for SRS/SBRT)',[],[], gd_pass};
else
    Title={'I. Geometric Distortion',  '(for non-SRS/SBRT)',[],[], gd_pass};
    limitlength = '0.2';
end
headerContent = {'Indicators','Measured Data (cm)','Reference Data (cm)','Abs.Diff. (cm)', 'Tolerance limit (cm)'};
bodyContent = {'A-B',num2str(dataAB,'%.2f'), num2str(refdataAB,'%.2f'), num2str(abs(refdataAB-dataAB),'%.2f'), limitlength;...
       'B-C',num2str(dataBC,'%.2f'), num2str(refdataBC,'%.2f'), num2str(abs(refdataBC-dataBC),'%.2f'), limitlength;...
       'C-A',num2str(dataCA,'%.2f'), num2str(refdataCA,'%.2f'), num2str(abs(refdataCA-dataCA),'%.2f'), limitlength;...
       'B-D',num2str(dataBD,'%.2f'), num2str(refdataBD,'%.2f'), num2str(abs(refdataBD-dataBD),'%.2f'), limitlength};

emptyContent = {[],[],[],[],[]};

ca = [Title; headerContent;bodyContent;emptyContent];
table1 = Table(ca);
table1.Style = {FontSize('13'),RowSep('single','black','1px')};
table1.Width = '100%';

if(strcmp(gd_pass,'Fail'))
    table1.Children(1,1).Children(1,5).Style = {HAlign('right'), Color("#FF0000")};
else
    table1.Children(1,1).Children(1,5).Style = {HAlign('right'), Color("#008000")};
end

for i=2:length(table1.Children)
    if(i==2)
        table1.Children(1,i).Style = {HAlign('center'), BackgroundColor("LightGray")};
    else
        table1.Children(1,i).Style = {HAlign('center')};
    end
end
add(rpt,table1);
add(rpt,text0);

% Print Uniformity Test
if(get(handles.radiobutton3,'Value'))
    un_pass = 'Pass';
else
    un_pass = 'Fail';
end
bodyContent = {'II. Uniformity and Noist Tests ', un_pass};
emptyContent = {[],[]};
ca = bodyContent;
table2 = Table(ca);
table2.Style = {FontSize('13')};
table2.Width = '100%';
if(strcmp(un_pass,'Fail'))
    table2.Children(1,1).Children(1,2).Style = {HAlign('right'), Color("#FF0000")};
else
    table2.Children(1,1).Children(1,2).Style = {HAlign('right'), Color("#008000")};
end
add(rpt,table2);


un_tol = get(handles.edit6,'String');
Title={'II.1. Uniformity Test', [], []};
headerContent = {'Max Measured HU Diff.','Max Ref. HU Diff','Tolerance limits.'};
bodyContent = {num2str(maxDiff_Uniform,'%.2f'), num2str(maxDiff_Uniform_ref,'%.2f'), un_tol};
emptyContent = {[],[],[]};



ca = [Title; headerContent;bodyContent;emptyContent];
table3 = Table(ca);
table3.Style = {FontSize('13'),RowSep('single','black','1px')};
table3.Width = '100%';
table3.Children(1,1).Children(1,3).Style = {HAlign('right')};

for i=2:length(table3.Children)
    if(i==2)
        table3.Children(1,i).Style = {HAlign('center'), BackgroundColor("LightGray")};
    else
        table3.Children(1,i).Style = {HAlign('center')};
    end
end
add(rpt,table3);

% Print Noise Test
noise_tol = get(handles.edit4,'String');

Title={'II.2. Noise Test', [], [], []};
headerContent = {[], 'Measured Noise (¥ò)','Ref. Noise (¥ò)','Tolerance limits (¥ò)'};
bodyContent1 = {'Small ROI',num2str(maxDiff_Noise,'%.2f'), num2str(maxDiff_Noise_Ref,'%.2f'), noise_tol};
bodyContent2 = {'Big ROI',num2str(maxDiff_Noise_BigCircle,'%.2f'), num2str(maxDiff_Noise_BigCircle_Ref,'%.2f'), noise_tol};
emptyContent = {[],[],[],[]};

ca = [Title; headerContent;bodyContent1;bodyContent2;emptyContent];
table4 = Table(ca);
table4.Style = {FontSize('13'), RowSep('single','black','1px')};
table4.Width = '100%';

for i=2:length(table4.Children)
    if(i==2)
        table4.Children(1,i).Style = {HAlign('center'), BackgroundColor("LightGray")};
    else
        table4.Children(1,i).Style = {HAlign('center')};
    end
end
add(rpt,table4);

add(rpt,text0);

% Print IVDT contrast Test
if(get(handles.radiobutton5,'Value'))
    IVDT_pass = 'Pass';
else
    IVDT_pass = 'Fail';
end
bodyContent = {'III. Constrast Test', IVDT_pass};
emptyContent = {[],[]};
ca = [bodyContent;emptyContent];
table5 = Table(ca);
table5.Style={FontSize('13')};
table5.Width = '100%';
if(strcmp(IVDT_pass,'Fail'))
    table5.Children(1,1).Children(1,2).Style = {HAlign('right'), Color("#FF0000")};
else
    table5.Children(1,1).Children(1,2).Style = {HAlign('right'), Color("#008000")};
end
add(rpt,table5);

add(rpt,text0);

% Print Spatial Resolution 
if(get(handles.radiobutton7,'Value'))
    sr_pass = 'Pass';
else
    sr_pass = 'Fail';
end
bodyContent = {'IV. Spatial Resolution Test', sr_pass};
emptyContent = {[],[]};
ca = [bodyContent;emptyContent];
table6 = Table(ca);
table6.Style={FontSize('13')};
table6.Width = '100%';

if(strcmp(sr_pass,'Fail'))
    table6.Children(1,1).Children(1,2).Style = {HAlign('right'), Color("#FF0000")};
else
    table6.Children(1,1).Children(1,2).Style = {HAlign('right'), Color("#008000")};
end
add(rpt,table6);


% Make Screenshot image
img = getframe(gcf);
imwrite(img.cdata, 'QA.png');
imageObj = Image('QA.png');
imageObj.Style = {ScaleToFit};

% Add Screenshot image to PDF
add(rpt,imageObj);

% Print PDF for QA report
rptview(rpt);
delete('QA.png');
set(handles.figure1, 'pointer', oldpointer);

end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function edit1_Callback(~, ~, ~)
end

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(~, ~, ~)
end

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(~, ~, ~)
end

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(~, ~, ~)
end
% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(~, ~, ~)
end

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(~, ~, ~)
end

% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(~, ~, ~)
end

% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(~, ~, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(handles.radiobutton7,'Value'))
    set(handles.radiobutton8,'Value',0);
    set(handles.radiobutton8,'FontWeight','normal');
    set(handles.radiobutton7,'FontWeight','bold');
else
    set(handles.radiobutton8,'Value',1);
    set(handles.radiobutton7,'FontWeight','normal');
    set(handles.radiobutton8,'FontWeight','bold');
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton7
end

% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(handles.radiobutton8,'Value'))
    set(handles.radiobutton7,'Value',0);
    set(handles.radiobutton7,'FontWeight','normal');
    set(handles.radiobutton8,'FontWeight','bold');
else
    set(handles.radiobutton7,'Value',1);
    set(handles.radiobutton8,'FontWeight','normal');
    set(handles.radiobutton7,'FontWeight','bold');
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton8
end



function windows_initial(handles)
global myDocsFolder;
% Axes
A=zeros(512,512);
imshow(A,'Parent', handles.axes1);
imshow(A,'Parent', handles.axes2);
imshow(A,'Parent', handles.axes3);
imshow(A,'Parent', handles.axes4);
imshow(A,'Parent', handles.axes5);
imshow(A,'Parent', handles.axes6);
imshow(A,'Parent', handles.axes7);

% Table1
materialData=zeros(20,4);
set(handles.uitable1, 'Data', materialData);

set(handles.radiobutton1,'Value',0);
set(handles.radiobutton1,'ForegroundColor',[0 128/255 0]);
set(handles.radiobutton1,'Enable','off');
set(handles.radiobutton2,'Value',1);
set(handles.radiobutton2,'ForegroundColor','r');
set(handles.radiobutton2,'Enable','off');
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton3,'ForegroundColor',[0 128/255 0]);
set(handles.radiobutton3,'Enable','off');
set(handles.radiobutton4,'Value',1);
set(handles.radiobutton4,'ForegroundColor','r');
set(handles.radiobutton4,'Enable','off');
set(handles.radiobutton5,'Value',0);
set(handles.radiobutton5,'ForegroundColor',[0 128/255 0]);
set(handles.radiobutton5,'Enable','off');
set(handles.radiobutton6,'Value',1);
set(handles.radiobutton6,'ForegroundColor','r');
set(handles.radiobutton6,'Enable','off');
set(handles.radiobutton7,'Value',0);
set(handles.radiobutton7,'ForegroundColor',[0 128/255 0]);
set(handles.radiobutton8,'Value',1);
set(handles.radiobutton8,'FontWeight','bold');
set(handles.radiobutton8,'ForegroundColor','r');
set(handles.radiobutton9,'Value',0);
set(handles.radiobutton10,'Value',1);
set(handles.radiobutton10,'FontWeight','bold');
set(handles.edit5,'Enable','off');
set(handles.checkbox1,'Value',1);
set(handles.checkbox2,'Value',1);
set(handles.checkbox3,'Value',1);
set(handles.checkbox4,'Value',1);
set(handles.checkbox5,'Value',1);
set(handles.checkbox6,'Value',1);
set(handles.checkbox7,'Value',1);
set(handles.checkbox8,'Value',1);
set(handles.checkbox9,'Value',1);
set(handles.checkbox10,'Value',1);
set(handles.checkbox11,'Value',1);
set(handles.checkbox12,'Value',1);
set(handles.checkbox13,'Value',1);
set(handles.checkbox14,'Value',1);
set(handles.checkbox15,'Value',1);
set(handles.checkbox16,'Value',1);
set(handles.checkbox17,'Value',1);
set(handles.checkbox18,'Value',1);
set(handles.checkbox19,'Value',1);
set(handles.checkbox20,'Value',1);

% table initialization
tabledata = zeros(20,4);
tabledataHTML = reshape(strtrim(cellstr(num2str(tabledata(:),'%.2f'))), size(tabledata));
tabledataHTML = strcat('<html><div style="width: 90; text-align: center;">', tabledataHTML,'</div></html>');
set(handles.uitable1, 'Data', tabledataHTML);

% set qa parameters from qaparam.dat
fid=fopen([myDocsFolder 'qaparam.dat'],'r');
smallR=fscanf(fid,'%f',[1 1]);
bigR=fscanf(fid,'%f',[1 1]);
maxHUdiff_noise=fscanf(fid,'%f',[1 1]);
useDC=fscanf(fid,'%f',[1 1]);
useSRS=fscanf(fid,'%f',[1 1]);
maxHUdiff_water=fscanf(fid,'%f',[1 1]);
maxHUdiff_others=fscanf(fid,'%f',[1 1]);
fiducialHU=fscanf(fid,'%f',[1 1]);

fclose(fid);

set(handles.edit2,'String',num2str(smallR,'%d'));
set(handles.edit3,'String',num2str(bigR,'%d'));    
set(handles.edit4,'String',num2str(maxHUdiff_noise,'%d'));
set(handles.edit5,'String',num2str(maxHUdiff_water,'%d'));
set(handles.edit6,'String',num2str(maxHUdiff_others,'%d'));
set(handles.edit7,'String',num2str(fiducialHU,'%d'));
if(useDC==1)
    set(handles.radiobutton9,'Value',1);
    set(handles.radiobutton10,'Value',0);
else
    set(handles.radiobutton9,'Value',0);
    set(handles.radiobutton10,'Value',1);
end

if(useSRS==1)
    set(handles.radiobutton11,'Value',1);
    set(handles.radiobutton12,'Value',0);
else
    set(handles.radiobutton11,'Value',0);
    set(handles.radiobutton12,'Value',1);
end



end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(~, ~, handles)
global myDocsFolder;
global refdataAB;
global refdataBC;
global refdataCA;
global refdataBD;
global refdataMean;
global refdataStdev;
global refdataMaterial;

global dataAB;
global dataBC;
global dataCA;
global dataBD;
global dataMean;
global dataStdev;
global materialData;
global step2_flag;
global acquisitionDate;
if(step2_flag)
    %d = char(datetime('today'));
    d = acquisitionDate;
    [savefile,savepath]=uiputfile({'*.mbl','MVCT Baseline File (*.mbl)'},'Save - Select a File',[myDocsFolder 'tomo_baseline_' d]);
    fid=fopen([savepath savefile], 'w');
    fprintf(fid, '%.2f\n', dataAB);
    fprintf(fid, '%.2f\n', dataBC);
    fprintf(fid, '%.2f\n', dataCA);
    fprintf(fid, '%.2f\n', dataBD);
    for i=1:length(dataMean)
        fprintf(fid, '%.2f %.2f\n', dataMean(i),dataStdev(i));
    end
    cnt = 0;
    for i=1:20
        eval(['val = get(handles.checkbox' num2str(i) ',''Value'');']);
        if(val)
            cnt=cnt+1;
        end
    end
    fprintf(fid, '%d\n', cnt);
    for i=1:length(materialData)
        eval(['val = get(handles.checkbox' num2str(i) ',''Value'');']);
        if(val)
            fprintf(fid, '%.2f\n', materialData(i,1));
        end
    end
    fprintf(fid, acquisitionDate);
    fclose(fid);
    refdataAB=dataAB;
    refdataBC=dataBC;
    refdataCA=dataCA;
    refdataBD=dataBD;
    refdataMean=dataMean;
    refdataStdev=dataStdev;
    refdataMaterial=materialData;
    pushbutton2_Callback(handles,handles,handles);
else
    msgbox('Please click "Step2" before.', 'Save Data Error!');
end


end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(~, ~, handles)
global myDocsFolder;
global refdataAB;
global refdataBC;
global refdataCA;
global refdataBD;
global refdataMean;
global refdataStdev;
global refdataMaterial;
global refacquisitionDate;
[loadfile,loadpath]=uigetfile({'*.mbl','MVCT Baseline File (*.mbl)'},'Load - Select a Baseline File',myDocsFolder);
fid=fopen([loadpath loadfile],'r');
refdataAB=fscanf(fid,'%f',[1 1]);
refdataBC=fscanf(fid,'%f',[1 1]);
refdataCA=fscanf(fid,'%f',[1 1]);
refdataBD=fscanf(fid,'%f',[1 1]);
data=fscanf(fid,'%f %f',[2 6]);
data=data';
refdataMean=data(1:end,1);
refdataStdev=data(1:end,2);
cnt=fscanf(fid,'%f',[1 1]);
refMat = zeros(20,1);
refMat(1:cnt) = fscanf(fid,'%f',[1 cnt]);    
refdataMaterial=refMat;
refacquisitionDate = fscanf(fid,'%s');    

fclose(fid);

if(refdataAB~=0) 
    set(handles.text16,'String','Reference data is loaded.');
else
    set(handles.text16,'String','Reference data is not loaded.');
end
end



% --- Executes on button press in checkbox1.
function checkbox1_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox3.
function checkbox3_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox6.
function checkbox6_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox7.
function checkbox7_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox8.
function checkbox8_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox9.
function checkbox9_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox11.
function checkbox11_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox12.
function checkbox12_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox13.
function checkbox13_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox14.
function checkbox14_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox15.
function checkbox15_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox16.
function checkbox16_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox17.
function checkbox17_Callback(~, ~, ~)
end

% --- Executes on button press in checkbox18.
function checkbox18_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox19.
function checkbox19_Callback(~, ~, ~)
end


% --- Executes on button press in checkbox20.
function checkbox20_Callback(~, ~, ~)
end


function edit2_Callback(~, ~, ~)
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
% --- Executes on button press in radiobutton7.
function radiobutton9_Callback(~, ~, handles)
if(get(handles.radiobutton9,'Value'))
    set(handles.radiobutton10,'Value',0);
    set(handles.radiobutton10,'FontWeight','normal');
    set(handles.radiobutton9,'FontWeight','bold');
    set(handles.edit5,'Enable','on');
else
    set(handles.radiobutton10,'Value',1);
    set(handles.radiobutton9,'FontWeight','normal');
    set(handles.radiobutton10,'FontWeight','bold');
    set(handles.edit5,'Enable','off');
end
end

% --- Executes on button press in radiobutton8.
function radiobutton10_Callback(~, ~, handles)
if(get(handles.radiobutton10,'Value'))
    set(handles.radiobutton9,'Value',0);
    set(handles.radiobutton9,'FontWeight','normal');
    set(handles.radiobutton10,'FontWeight','bold');    
    set(handles.edit5,'Enable','off');
else
    set(handles.radiobutton9,'Value',1);
    set(handles.radiobutton10,'FontWeight','normal');
    set(handles.radiobutton9,'FontWeight','bold');
    set(handles.edit5,'Enable','on');
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton8
end


function edit4_Callback(~, ~, ~)
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit3_Callback(~, ~, ~)
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit5_Callback(~, ~, ~)
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit6_Callback(~, ~, ~)
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(~, ~, handles)
global myDocsFolder;
    % save new QA parameters to qaparam.dat file
    fid=fopen([myDocsFolder 'qaparam.dat'], 'w');
    fprintf(fid, '%d\n', uint16(str2double(get(handles.edit2,'String'))));
    fprintf(fid, '%d\n', uint16(str2double(get(handles.edit3,'String'))));
    fprintf(fid, '%.1f\n', str2double(get(handles.edit4,'String')));
    if(get(handles.radiobutton9,'Value'))
        fprintf(fid, '1\n');        
    else
        fprintf(fid, '0\n');
    end    
    if(get(handles.radiobutton11,'Value'))
        fprintf(fid, '1\n');        
    else
        fprintf(fid, '0\n');
    end    
    fprintf(fid, '%d\n',  uint16(str2double(get(handles.edit5,'String'))));    
    fprintf(fid, '%d\n',  uint16(str2double(get(handles.edit6,'String'))));
    fprintf(fid, '%d\n',  uint16(str2double(get(handles.edit7,'String'))));
    fclose(fid);
end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(~, ~, handles)
    global myDocsFolder;
    % save default QA parameters to qaparam.dat file   
    fid=fopen([myDocsFolder 'qaparam.dat'], 'w');    
    fprintf(fid, '5\n120\n70\n0\n1\n30\n50\n1400\n');    
    fclose(fid);
    % set default QA parameters
    set(handles.edit2,'String','5');
    set(handles.edit3,'String','120');
    set(handles.edit4,'String','70');    
    set(handles.edit5,'String','30');
    set(handles.edit6,'String','50');
    set(handles.edit7,'String','1400');
    set(handles.radiobutton9,'Value',0);    
    set(handles.radiobutton10,'Value',1);    
    set(handles.radiobutton11,'Value',1);
    set(handles.radiobutton12,'Value',0);    
end


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
if(get(handles.radiobutton11,'Value'))
    set(handles.radiobutton12,'Value',0);
    set(handles.radiobutton12,'FontWeight','normal');
    set(handles.radiobutton11,'FontWeight','bold');
else
    set(handles.radiobutton12,'Value',1);
    set(handles.radiobutton11,'FontWeight','normal');
    set(handles.radiobutton12,'FontWeight','bold');
end
end

% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
if(get(handles.radiobutton12,'Value'))
    set(handles.radiobutton11,'Value',0);
    set(handles.radiobutton11,'FontWeight','normal');
    set(handles.radiobutton12,'FontWeight','bold');
else
    set(handles.radiobutton11,'Value',1);
    set(handles.radiobutton12,'FontWeight','normal');
    set(handles.radiobutton11,'FontWeight','bold');
end
end



function edit7_Callback(~,~,~)
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
