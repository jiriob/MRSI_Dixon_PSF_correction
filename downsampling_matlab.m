clear all;

directory = '/Volumes/Home/zgung/Desktop/Current/Schlapata/';
field = 3;
no_cor = 0;
jmr = 0;
shft_ud = 0;
shft_lr = 0;
SNR_filter = 0.2; % it is 0.1 for phantom and 0.2 for in vivo
dm1 = 128;
dm2 = 128;
dm3 = 64;



% minarikova.lenka@gmail.com

% !!!!!!!!!!!!!!!!!! readme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% for working you need my other function called read_ascconv_lenk.m 
%   for reading parameters from the dicom
%   and you need a txt file with a set of SNR to evaluate, made with 
%   my other function called SNR_lenk.m
% directory = '~/Patient_name' - where is a directory called "Spec" with 
%   a dicom 3D CSI file
% field = 3 - e
% w_nois water pixel intensity noise threshold
% f_nois = fat noise threshold (ex. I know fat intensity is only above 50)
% s_nois: set to 1 if you want saved noise information to be loaded 
%   instead of again computing new data sets 0
% jmr = 1 if the data are from jmrui (usualy for phantom measurements)
% shft_ud: for shifted CSI: up -0.% down +0.%
% shft_lr: for shifted CSI: left -0.% right +0.%
% SNR_filter: change this value to 0.1 if you want 10% of SNR to be filter (for 
%    phantoms) and to 0.2 (for in vivo data)

% the output is maximal, mean value of all SNRs of Cho and a table 
%   with all SNRs in one row, all saved in txt files in Spec directory

%if there is this error: Undefined function or variable "w_ind".
% just remove the noise!

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% search for variables in spectroscopy file:
tic;
%directory_main = ('/Volumes/Home/zgung/Desktop/cholin/Kuktics/');
%cd(directory_main);
% List all data folders
%subs = dir(pwd);
% Remove non-folders, . and .. items from data folder list
%for k = length(subs):-1:1
    % extract subjects names
%    if ~subs(k).isdir
%        subs(k) = [];
%        continue
%    end
%    fname = subs(k).name;
%    if fname(1) == '.'
%        subs(k) = [];
%    end
%end
%%%
%subjects = subs;
%nextSub = 1;% <<<<<<<<< Here you can choose the data folder 
%directory = strcat(directory_main,subjects(nextSub).name,'/');
directory_main = directory;
directory_spec = strcat(directory,'Spec/');
cd(directory_spec);
disp(strcat('Processing:',directory));
%% <<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>
press_big = 1; % choose if you want more voxel (even that no totally inside the press box)
%cntrlfor1time = 1;
%field = 3; % 3 or 7 for different dicoms from different coils
%w_nois = 75 + 3 * 6.38; %read from osirix (mean value of noise + 3*STDEV)
%f_nois = 30 + 3 * 7.7;
%csishift = 1; %if the matrix is shifted of 25% (4.2) up

% read from spectroscopic file from the text part
% beginnig by "### ASCCONV BEGIN ###" and ending by "### ASCCONV END ###"
slices = dir(pwd);
for k = length(slices):-1:1 % find .IMA file
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    elseif strcat(fname(end-2),fname(end-1),fname(end)) ~= strcat('IMA')
        slices(k) = [];
    end
end
voxel = read_ascconv_lenk(slices(1).name); % read parrameter from spectroscopy dicom

% check coil elements:
if field == 3 % anyway we are not measuring dixons on 7 T
    if voxel.coilel3 ~= 0
        disp('!!Both coils ON!!');
    else
        if strcat(voxel.coilel1(4)) == strcat('R') && ...
                strcat(voxel.coilel2(4)) == strcat('R')
            disp('Right coil elements: ON');
        end
        if strcat(voxel.coilel1(4)) == strcat('L') && ...
                strcat(voxel.coilel2(4)) == strcat('L')
            disp('Left coil elements: ON');
        end
    end  
end
%% move to dixom folder: water first
disp('Processing PSF');
cd(strcat(directory,'Dixon_1.0iso_PAT2_v2_W/'));
slices = dir(pwd);
for k = length(slices):-1:1
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    end
end
%% determine CSI-in-press parameters:
voxel.size_z = voxel.FoV_z / voxel.number_z; % the size of 1 voxel after zero filling in mm
voxel.size_y = voxel.FoV_y / voxel.number_y;
voxel.size_x = voxel.FoV_x / voxel.number_x;
% number of steps in pressbox after zero filling in each direction / 2:
if press_big == 1
    voxel.step_x = voxel.number_x / 2 - ceil((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x); %include also voxels touching the pressbox fringes
    voxel.step_y = voxel.number_y / 2 - ceil((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y);
    voxel.step_z = voxel.number_z / 2 - ceil((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z);
elseif press_big == 0
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x) - 1; %include only voxels inside the PRESS box
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y) - 1;
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z) - 1;
else
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x); %include voxels with PRESSbox fringes inside
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y);
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z);
end

nfo1 = dicominfo(slices(1).name);
%% Import the images, center them and regrid
% make the center of coordinates the centre of image: horizontaly
n = numel(slices);
imgs = cell(1,1);
for ii = 1:n
    imgs{1,1}(:,:,ii) = double(dicomread(slices(ii).name));
end
%% simulate shorter measurement: cut the center of k-space to half
% add zeros so you have the same resolution
fft_imgs{1,1} = fftshift(fftn(imgs{1,1}));
fft_imgs{2,1} = fft_imgs{1,1}((nfo1.Rows - dm1) / 2:(nfo1.Rows + dm1) / 2 - 1,(nfo1.Columns - dm2) / 2:(nfo1.Columns + dm2) / 2 - 1,...
    (n - dm3) / 2:(n + dm3) / 2 - 1);
fft_imgs{3,1} = abs(ifftn(padarray(fft_imgs{2,1},[double(nfo1.Rows - dm1) / 2,double(nfo1.Columns - dm2) / 2,double(n - dm3) / 2])));
%
[XI,YI,ZI] = meshgrid(1:(floor(10000 * 1 / nfo1.PixelSpacing(1,1)) / 10000):length(imgs{1,1}(1,:,1)),...
        1:(floor(10000 * 1 / nfo1.PixelSpacing(2,1)) / 10000):length(imgs{1,1}(:,1,1)),...
        1:(floor(10000 * 1 / nfo1.SliceThickness) / 10000):length(imgs{1,1}(1,1,:)));
    imgs{2,1} = ba_interp3(fft_imgs{3,1},XI,YI,ZI,'nearest'); % fast interpolation to 1x1x1 mm3
    [YI,XI,~] = size(imgs{2,1});

%% make the center of coordinates the centre of each axial slice: horizontaly
if -nfo1.ImagePositionPatient(1,1) > XI -  ...
        -nfo1.ImagePositionPatient(1,1) == 1
        % add zeros to the end
        imgs{2,1}(:,XI + 1:-round(2 * ...
            nfo1.ImagePositionPatient(1,1)),:) = 0;
else
        % move image right and add zeros to the begging
        imgs{2,1}(:,XI - -round(2 * nfo1.ImagePositionPatient(1,1)) + ...
            1:2 * XI - -round(2 * nfo1.ImagePositionPatient(1,1)),:) = ...
            imgs{2,1}(:,1:XI,:);
        imgs{2,1}(:,1:XI - -round(2 * nfo1.ImagePositionPatient(1,1)),:) = 0;
end
% now, center image in vertical dimension
if -nfo1.ImagePositionPatient(2,1) > YI - -nfo1.ImagePositionPatient(2,1) == 1
        %add zeros to the end
        imgs{2,1}(YI + 1:-round(2 * ...
            nfo1.ImagePositionPatient(2,1)),:,:) = 0;
else
        % move image down and add zeros to the begging
        imgs{2,1}(YI - -round(2 * nfo1.ImagePositionPatient(2,1)) ...
            + 1:2 * YI - -round(2 * nfo1.ImagePositionPatient(2,1)),:,:) = ...
            imgs{2,1}(1:YI,:,:);
        imgs{2,1}(1:YI - -round(2 * nfo1.ImagePositionPatient(2,1)),:,:) = 0;
end
%% other csi parameters:
% determine the first slice and last slices = determine the press box transversal beginning:
voxel.FoV_z_1 = round(voxel.fov_cntr_z + voxel.p_fov_z / 2); % the real position in number of pixels in the image
voxel.FoV_z_2 = round(voxel.fov_cntr_z - voxel.p_fov_z / 2);
voxel.csi_slice_1 = round((nfo1.ImagePositionPatient(3,1) - voxel.FoV_z_1) + 1) + 2; %1st - fhead
voxel.csi_slice_last = round((nfo1.ImagePositionPatient(3,1) - voxel.FoV_z_2)) + 2; %2nd - from feet
imgs{3,1} = imgs{2,1}(:,:,voxel.csi_slice_1:voxel.csi_slice_last); % remove the additional slices    

% look if the csi is rotated about an angle:
voxel.angle_deg = radtodeg(voxel.angle);
%%
% the center point of csi in press box in rotated images has coordinates from the
% left upper corner:
if voxel.angle > 0.0 % clockwise
    voxel.fov_cntr_x = voxel.fov_cntr_x + 4;
    voxel.fov_cntr_y = voxel.fov_cntr_y + 2; % plus posunie ratio hore
elseif voxel.angle < -0.0 % anticlockwise
    voxel.fov_cntr_x = voxel.fov_cntr_x + 4;
    voxel.fov_cntr_y = voxel.fov_cntr_y + 3; % plus posunie ratio hore
else
    voxel.fov_cntr_x = voxel.fov_cntr_x + 4;
    voxel.fov_cntr_y = voxel.fov_cntr_y + 3; % plus posunie ratio hore
end
voxel.fov_cntr_x_rttd = numel(imgs{3,1}(1,:,1)) / 2 - (-voxel.fov_cntr_x) * cos(-voxel.angle) - ...
     (-voxel.fov_cntr_y) * sin(-voxel.angle) + 0;
voxel.fov_cntr_y_rttd = numel(imgs{3,1}(:,1,1)) / 2 - (-voxel.fov_cntr_y) * cos(-voxel.angle) + ...
     (-voxel.fov_cntr_x) * sin(-voxel.angle) + 0;

% and width and hight:
voxel.fov_x1 = floor(voxel.fov_cntr_x_rttd - voxel.p_fov_x / 2);
voxel.fov_y1 = floor(voxel.fov_cntr_y_rttd - voxel.p_fov_y / 2);
% read image and turn it, cut it
for ii = 1:numel(imgs{3,1}(1,1,:))
   imgs{4,1}(:,:,ii) = imcrop(imrotate(imgs{3,1}(:,:,ii),-voxel.angle_deg,...
       'nearest','crop'),[voxel.fov_x1 voxel.fov_y1 (voxel.p_fov_x - 1) (voxel.p_fov_y - 1)]);
% rotate default pictures:
   imgs{3,1}(:,:,ii) = imrotate(imgs{3,1}(:,:,ii),-voxel.angle_deg,...
       'nearest','crop');
end
%% tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
%%%%%%% move to dixom folder: fat
cd(strcat(directory,'Dixon_1.0iso_PAT2_v2_F/'));
slices = dir(pwd);
for k = length(slices):-1:1
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    end
end
% make the center of coordinates the centre of image: horizontaly
n = numel(slices);
for ii = 1:n
    imgs{1,2}(:,:,ii) = double(dicomread(slices(ii).name));
end
fft_imgs{1,2} = fftshift(fftn(imgs{1,2}));
fft_imgs{2,2} = fft_imgs{1,2}((nfo1.Rows - dm1) / 2:(nfo1.Rows + dm1) / 2 - 1,(nfo1.Columns - dm2) / 2:(nfo1.Columns + dm2) / 2 - 1,...
    (n - dm3) / 2:(n + dm3) / 2 - 1);
fft_imgs{3,2} = abs(ifftn(padarray(fft_imgs{2,2},[double(nfo1.Rows - dm1) / 2,double(nfo1.Columns - dm2) / 2,double(n - dm3) / 2])));

[XI,YI,ZI] = meshgrid(1:1 / nfo1.PixelSpacing(1,1):length(imgs{1,2}(1,:,1)),...
    1:1 / nfo1.PixelSpacing(2,1):length(imgs{1,2}(:,1,1)),...
    1:1 / nfo1.SliceThickness:length(imgs{1,2}(1,1,:)));
imgs{2,2} = ba_interp3(fft_imgs{3,2},XI,YI,ZI,'nearest'); % fast interpolation to 1x1x1 mm3
[YI,XI,~] = size(imgs{2,2});
% make the center of coordinates the centre of image: horizontaly
if -nfo1.ImagePositionPatient(1,1) > XI -  ...
        -nfo1.ImagePositionPatient(1,1) == 1
        % add zeros to the end
        imgs{2,2}(:,XI + 1:-round(2 * ...
            nfo1.ImagePositionPatient(1,1)),:) = 0;
else
        % move image right and add zeros to the begging
        imgs{2,2}(:,XI - -round(2 * nfo1.ImagePositionPatient(1,1)) + ...
            1:2 * XI - -round(2 * nfo1.ImagePositionPatient(1,1)),:) = ...
            imgs{2,2}(:,1:XI,:);
        imgs{2,2}(:,1:XI - -round(2 * nfo1.ImagePositionPatient(1,1)),:) = 0;
end
% now, center image in vertical dimension
if -nfo1.ImagePositionPatient(2,1) > YI - -nfo1.ImagePositionPatient(2,1) == 1
        %add zeros to the end
        imgs{2,2}(YI + 1:-round(2 * ...
            nfo1.ImagePositionPatient(2,1)),:,:) = 0;
else
        % move image right and add zeros to the begging
        imgs{2,2}(YI - -round(2 * nfo1.ImagePositionPatient(2,1)) ...
            + 1:2 * YI - -round(2 * nfo1.ImagePositionPatient(2,1)),:,:) = ...
            imgs{2,2}(1:YI,:,:);
        imgs{2,2}(1:YI - -round(2 * nfo1.ImagePositionPatient(2,1)),:,:) = 0;
end
imgs{3,2} = imgs{2,2}(:,:,voxel.csi_slice_1:voxel.csi_slice_last); % remove the additional slices
% read image and turn it, cut it
for ii = 1:numel(imgs{3,2}(1,1,:))
   imgs{4,2}(:,:,ii) = imcrop(imrotate(imgs{3,2}(:,:,ii),-voxel.angle_deg,...
       'nearest','crop'),[voxel.fov_x1 voxel.fov_y1 (voxel.p_fov_x - 1) (voxel.p_fov_y - 1)]);
% rotate default pictures:
   imgs{3,2}(:,:,ii) = imrotate(imgs{3,2}(:,:,ii),-voxel.angle_deg,...
       'nearest','crop');
end

%% %%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(strcat(directory,'Spec/'));
ova(1,1) = numel(imgs{4,1}(:,1,1)); % x
ova(2,1) = numel(imgs{4,1}(1,:,1)); % y
ova(3,1) = numel(imgs{4,1}(1,1,:)); % z

threed(:,1) = double(reshape(imgs{4,1},[],1));
threed(:,2) = double(reshape(imgs{4,2},[],1));
for i = 1:ova(3,1) %do I need this???
    threed(((i-1) * ova(1,1) * ova(2,1)) + 1: i * ova(1,1) * ova(2,1),3) = i;
end
disp('Segmenting...');
%% segmentation alone
data = threed(:,1:2);
% how to distinquish water cluster from noise and fat using 2d histogram
[hist,biny] = hist3(data,[128,128]);
biny_a(:,1) = biny{1,1};
biny_a(:,2) = biny{1,2};
clear biny;
iii = 100;
figure(1);
subplot(221);
imagesc(hist);
caxis([0 100]);
% filter noise in the histogram
noiseFilter = imhmin(hist,2);
%# Gaussian-filter the image:
gaussFilter = fspecial('gaussian',[64 64],5);  %# Create the filter
filteredData = imfilter(noiseFilter,gaussFilter);
subplot(222);
imagesc(filteredData);
title('Gaussian-filtered image');
caxis([0 100]);
% # Perform a morphological close operation:
closeElement = strel('disk',4);  %# Create a disk-shaped structuring element
closedData = imclose(filteredData,closeElement);
subplot(223);
imagesc(closedData);
title('Closed image');
caxis([0 100]);
% find centers
mxma = imregionalmax(closedData);
subplot(224);
imagesc(mxma);
[X,Y] = find(mxma==max(mxma(:))); 
% assign bins to the max points
[w.X,w.ind] = max(X);
w.Y = Y(w.ind);
[f.Y,f.ind] = max(Y);
f.X = X(f.ind);
% find the "minimum" between the two maxima (or even three)
syms a k
[sol_a, sol_k] = solve(k * w.X + a == w.Y, k * f.X + a == f.Y);
yy = ones(size((f.X + 1):(w.X - 1)));
dif = yy;
i = 0;
if closedData(f.X,f.Y) > closedData(w.X,w.Y)
    for ii = (f.X + 1):(w.X - 1) % pouzi indexy radsej
        i = i + 1;
        yy(i) = round(double(sol_k) * ii + double(sol_a));
        % a teraz pekne sprav rozdiely
        dif(i) = closedData(f.X,f.Y) - closedData(ii,yy(i));
    end
else
    disp('Water peak is higher than fat?')
end
ii = (f.X + 1):(w.X - 1);
[~,coor_min.abs] = max(dif);
coor_min.X = ii(coor_min.abs);
coor_min.Y = yy(coor_min.abs);
% bins for the water center and for "minimum"
%f.wii = biny_a(f.X,1); % water images intensities
%f.fii = biny_a(f.Y,2); % fat image intensities
w.wii = biny_a(w.X,1);
w.fii = biny_a(w.Y,2);
coor_min.wii = biny_a(coor_min.X,1);
coor_min.fii = biny_a(coor_min.Y,2);
% solve a equation: a circle with the center in water maximum and radius of the "minimum"
% (x - w.wii)^2 + (y - w.fii)^2 = sqrt((coor_min.wii - w.wii)^2 + (coor_min.fii - w.fii)^2))
% second equation: y = x (where x is the "minimum")
for i = 1:length(threed)
    if threed(i,1) < w.wii % for voxels with smaller intensity on water images than coor_min.wii
        if (threed(i,1) - w.wii)^2 + (threed(i,2) - w.fii)^2 > (coor_min.wii - w.wii)^2 + (coor_min.fii - w.fii)^2
            threed(i,4) = 0;
            threed(i,5) = 0;
        else
            threed(i,4) = threed(i,1);
            threed(i,5) = 1;
        end
    else
        if  threed(i,2) > threed(i,1) - ...
                (w.wii - (sqrt((coor_min.wii - w.wii)^2 + (coor_min.fii - w.fii)^2) + w.fii))
            threed(i,4) = 0;
            threed(i,5) = 0;
        else
            threed(i,4) = threed(i,1);
            threed(i,5) = 1;
        end
    end
end
% gscatter(threed(:,1),threed(:,2),threed(:,4),'br','xo'), axis tight
figure(2);
gscatter(threed(:,1),threed(:,2),threed(:,5),'rg','.','','off'), axis tight
%% define which cluster is water and fat
w_ind = 1; f_ind = 0;
%% make maps with (real) 1 and 0 values
% make a segmented image: 
sgmnts{1,1} = reshape(threed(:,4),ova(1,1),ova(2,1),ova(3,1));
%sgmnts{1,2} = reshape(threed_(:,5),ova(1,1),ova(2,1),ova(3,1)); %segmented images for fat as 1
% for shifted CSI, you need to shift the ratio information:
if shft_ud < 0 || shft_ud > 0
    voxel.img_s = size(sgmnts{1,1});
    shft_ud = round(shft_ud * 100) / 100;
    [XI,YI,ZI] = meshgrid(1:voxel.img_s(2),...
    1:(voxel.img_s(1) - 0.95) / (voxel.img_s(1) * 10):voxel.img_s(1),...
    1:voxel.img_s(3));
    sgmnts{1,3} = ba_interp3(sgmnts{1,1},XI,YI,ZI,'linear');
%    sgmnts{1,4} = ba_interp3(sgmnts{1,2},XI,YI,ZI,'linear');
    % the shift:
    % in Y direction:
    [YI,~,~] = size(sgmnts{1,3});
    knst = round(10 * abs(shft_ud) * (voxel.FoV_y / voxel.notinterpfov_y));
    if shft_ud < 0 % shift it up
        sgmnts{1,3}(knst + 1:YI,:,:) = sgmnts{1,3}(1:YI - knst,:,:);
        sgmnts{1,3}(1:knst,:,:) = 0;
%        sgmnts{1,4}(knst + 1:YI,:,:) = sgmnts{1,4}(1:YI - knst,:,:);
%        sgmnts{1,4}(1:knst,:,:) = 0;
    end
    if shft_ud > 0 % shift it down
        sgmnts{1,3}(YI + 1:YI + knst,:,:) = 0;
        sgmnts{1,3}(1:YI,:,:) = sgmnts{1,3}(knst + 1:YI + knst,:,:);
%        sgmnts{1,4}(YI + 1:YI + knst,:,:) = 0;
%        sgmnts{1,4}(1:YI,:,:) = sgmnts{1,4}(knst + 1:YI + knst,:,:);
    end
    % interpolate back:
    [XI,YI,ZI] = meshgrid(1:voxel.img_s(2),1:voxel.img_s(1),1:voxel.img_s(3));
    [X,Y,Z] =  meshgrid(1:voxel.img_s(2),...
    1:(voxel.img_s(1) - 0.95) / (voxel.img_s(1) * 10):voxel.img_s(1),...
    1:voxel.img_s(3));
    sgmnts{1,1} = ba_interp3(X,Y,Z,sgmnts{1,3},XI,YI,ZI,'linear');
%    sgmnts{1,2} = ba_interp3(X,Y,Z,sgmnts{1,4},XI,YI,ZI,'linear');
    clear sgmnts{1,3};
%    clear sgmnts{1,4};
end
if shft_lr < 0 || shft_lr > 0
    voxel.img_s = size(sgmnts{1,1});
    shft_lr = round(shft_lr * 100) / 100;
    % make it 10 times bigger in X direction
    [XI,YI,ZI] = meshgrid(1:(voxel.img_s(2) - ...
    0.95) / (voxel.img_s(2) * 10):voxel.img_s(2),1:voxel.img_s(1),...
    1:voxel.img_s(3));
    sgmnts{1,3} = ba_interp3(sgmnts{1,1},XI,YI,ZI,'linear');
%    sgmnts{1,4} = ba_interp3(sgmnts{1,2},XI,YI,ZI,'linear');
    % the shift:
    % in X direction:
    [~,XI,~] = size(sgmnts{1,3});
    knst = round(10 * abs(shft_lr) * (voxel.FoV_x / voxel.notinterpfov_x));
    if shft_lr < 0 % shift it left
        sgmnts{1,3}(:,knst + 1:XI,:) = sgmnts{1,3}(:,1:XI - knst,:);
        sgmnts{1,3}(:,1:knst,:) = 0;
%        sgmnts{1,4}(:,knst + 1:XI,:) = sgmnts{1,4}(:,1:XI - knst,:);
%        sgmnts{1,4}(:,1:knst,:) = 0;
    end
    if shft_lr > 0 % shift it right
        sgmnts{1,3}(:,XI + 1:XI + knst,:) = 0;
        sgmnts{1,3}(:,1:XI,:) = sgmnts{1,3}(:,knst + 1:XI + knst,:);
%        sgmnts{1,4}(:,XI + 1:XI + knst,:) = 0;
%        sgmnts{1,4}(:,1:XI,:) = sgmnts{1,4}(:,knst + 1:XI + knst,:);
    end
    % interpolate back:
    [XI,YI,ZI] = meshgrid(1:voxel.img_s(2),1:voxel.img_s(1),1:voxel.img_s(3));
    [X,Y,Z] = meshgrid(1:(voxel.img_s(2) - ...
    0.95) / (voxel.img_s(2) * 10):voxel.img_s(2),1:voxel.img_s(1),...
    1:voxel.img_s(3));
    sgmnts{1,1} = ba_interp3(X,Y,Z,sgmnts{1,3},XI,YI,ZI,'linear');
%    sgmnts{1,2} = ba_interp3(X,Y,Z,sgmnts{1,4},XI,YI,ZI,'linear');
    clear sgmnts{1,3};
%    clear sgmnts{1,4};
end
%% add zeros to the map, so it's a cube with dimensions of 12x12x12 MRS voxels:
sgmnts{2,1} = zeros(voxel.rm,voxel.rm,voxel.rm); % make a 120x120x120 matrix composed of zeros
%sgmnts{2,2} = zeros(voxel.rm,voxel.rm,voxel.rm);
sgmnts{2,1}(round((voxel.rm / 2) + 1 - ova(1,1) / 2):round((voxel.rm / 2) + ...
    ova(1,1) / 2),round((voxel.rm / 2) + 1 - ova(2,1) / 2):round((voxel.rm / 2) + ...
    ova(2,1) / 2),round((voxel.rm / 2) + 1 - ova(3,1) / 2):round((voxel.rm / 2) + ...
    ova(3,1) / 2)) = sgmnts{1,1}(:,:,:);
%sgmnts{2,2}(round((voxel.rm / 2) + 1 - ova(1,1) / 2):round((voxel.rm / 2) + ...
%    ova(1,1) / 2),round((voxel.rm / 2) + 1 - ova(2,1) / 2):round((voxel.rm / 2) + ...
%    ova(2,1) / 2),round((voxel.rm / 2) + 1 - ova(3,1) / 2):round((voxel.rm / 2) + ...
%    ova(3,1) / 2)) = sgmnts{1,2}(:,:,:);
sgmnts{3,1} = fftshift(fftn(sgmnts{2,1})); % fft of the map
%sgmnts{3,2} = fftshift(fftn(sgmnts{2,2}));
% imagesc(sgmnts{2,1}(:,:,61)); axis equal
% %%
% %figure('units','normalized','position',[.05 .05 .6 .935]);
% imagesc(log(abs(sgmnts{3,1}(:,:,1)))); % to see the fft image of 1 slice
sgmnts_w = sgmnts{2,1};
save(strcat(nfo1.PatientName.FamilyName,'_sgmnts_w.mat'),'sgmnts_w');
clear sgmnts_save;

%% elliptical k-space !!!
sgmnts{4,1} = sgmnts{3,1}(((voxel.rm / 2) - 5):((voxel.rm / 2) + ...
    6),((voxel.rm / 2) - 5):((voxel.rm / 2) + 6),((voxel.rm / 2) - 5):((voxel.rm / 2) + 6));
%sgmnts{4,2} = sgmnts{3,2}(((voxel.rm / 2) - 5):((voxel.rm / 2) + ...
%    6),((voxel.rm / 2) - 5):((voxel.rm / 2) + 6),((voxel.rm / 2) - 5):((voxel.rm / 2) + 6));
%%
sgmnts{5,1} = sgmnts{4,1};
%sgmnts{5,2} = sgmnts{4,2};
pts_max = 12;
avg_max = 2;
if (mod(pts_max,2))
   kpts = pts_max/2;
else
   kpts = pts_max/2-1;
end
%elliptical_filter = zeros(pts_max,pts_max,pts_max);
for avg=1:avg_max
   for i=1:pts_max
       for j=1:pts_max
           for k=1:pts_max
%               dist=sqrt(((i-pts32/2)*(i-pts32/2)+(j-pts32/2)*(j-pts32/2)+(k-pts32/2)*(k-pts32/2)))/(pts32/2);
               dist=sqrt(((i-pts_max/2)*(i-pts_max/2)+(j-pts_max/2)*(j-pts_max/2)+(k-pts_max/2)*(k-pts_max/2)))/(kpts);
               if ( dist <= 1 )
                   sgmnts{5,1}(i,j,k) = (floor(0.5+(avg-1)*(0.5+0.5*cos(pi*dist))+1)) * sgmnts{5,1}(i,j,k);
%                   sgmnts{5,2}(i,j,k) = (floor(0.5+(avg-1)*(0.5+0.5*cos(pi*dist))+1)) * sgmnts{5,2}(i,j,k);
               else 
                   sgmnts{5,1}(i,j,k) = 0;
%                   sgmnts{5,2}(i,j,k) = 0;
               end
           end
       end
   end
end
%% hamming filter:
% generate 1D filter with CSI matrix resolution:
w1 = hamming(11);
[x,y,z] = meshgrid(-5:1:5);
r = sqrt(x.^2 + y.^2 + z.^2);
w = zeros(size(r));
w(r<=5) = interp1(linspace(-5,5,11),w1,r(r<=5));
%imagesc(w(:,:,6));
% fill the matrix to 12x12x12
w = padarray(w,[1 1 1],'post');
sgmnts{6,1} = sgmnts{5,1} .* w; % multiply the hamming filter with fourier transform
%sgmnts{6,2} = sgmnts{5,2} .* w;
%% zero filling to 16x16x16
sgmnts{7,1} = abs(ifftn(padarray(sgmnts{6,1},[2 2 2])));
%sgmnts{7,2} = abs(ifftn(padarray(sgmnts{6,2},[2 2 2])));
%% scale it!
ttl = max(max(max(sgmnts{7,1}))); % maximum value, not sure if this is the brightess idea
l = voxel.number_x / 2 - voxel.step_x;
m = voxel.number_y / 2 - voxel.step_y;
n = voxel.number_z / 2 - voxel.step_z;
for k = 1:voxel.step_z * 2
    for i = 1:voxel.step_x * 2
        for j = 1:voxel.step_y * 2
            sgmnts{10,1}(j,i,k) = sgmnts{7,1}(m + j,l + i,n + k) / ttl;
        end
    end
end

rtio = sgmnts{10,1}; % ./ (sgmnts{10,1} + sgmnts{10,2});

%% %%%%%%%%%%%%%%%% import the results from text file: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(directory_spec);
% import from simple text file, only amplitude data from jmrui text result
% file named Output_choSNR.txt
fid = fopen('Output_choSNR.txt','r');
mtrx1 = textscan(fid,'%f'); % this is right? not just '%f' ??
fclose(fid);
%% you need to cut the images in voxels dimensions (not in press box!)
c.cut_width = round(2 * voxel.step_y * voxel.size_y); % width in pixels of the cuted image
c.cut_height = round(2 * voxel.step_x * voxel.size_x); % height -||-
c.cut_depth = round(2 * voxel.step_z * voxel.size_z); % depth -||-
c.pix_width = voxel.step_y * 2; % number of voxels in the pressbox - x axis
c.pix_height = voxel.step_x * 2; % -||- - y axis
c.pix_depth = voxel.step_z * 2; % -||- - z axis
voxel.FoV_z_1 = round(voxel.fov_cntr_z + c.cut_depth / 2); % the real position in number of pixels in the image
voxel.FoV_z_2 = round(voxel.fov_cntr_z - c.cut_depth / 2);
voxel.csi_slice_1 = round((nfo1.ImagePositionPatient(3,1) - voxel.FoV_z_1) + 1) + 0; %1st - fhead
voxel.csi_slice_last = round((nfo1.ImagePositionPatient(3,1) - voxel.FoV_z_2)) + 0; %2nd - from feet
voxel.fov_x1 = floor(voxel.fov_cntr_x_rttd - voxel.step_x * voxel.size_x + round(shft_lr * (voxel.FoV_x / voxel.notinterpfov_x))); % voxel values are beeing saved for psf_pics.m!
voxel.fov_y1 = floor(voxel.fov_cntr_y_rttd - voxel.step_y * voxel.size_y + round(shft_ud * (voxel.FoV_y / voxel.notinterpfov_y)));

% cut the voxels-images only in Z direction:
imgs{5,1} = imgs{2,1}(:,:,voxel.csi_slice_1:voxel.csi_slice_last); % remove the additional slices
for ii = 1:numel(imgs{5,1}(1,1,:))
% rotate default pictures:
   imgs{5,1}(:,:,ii) = imrotate(imgs{5,1}(:,:,ii),-voxel.angle_deg,...
       'nearest','crop');
end
imgs_w = imgs{5,1};
save(strcat(nfo1.PatientName.FamilyName,'_imgs_w.mat'),'imgs_w');
clear imgs_w;

aa = 0; % Z direction in siemens MRSI is opposite to image Z direction
if jmr == 1
    for i = 1:voxel.step_x * 2
        for j = 1:voxel.step_y * 2
            for k = voxel.step_z * 2:-1:1
                aa = aa + 1;
                mtrx2(j,i,k) = mtrx1{1,1}(aa,1);
            end
        end
    end
else
    for k = voxel.step_z * 2:-1:1
        for j = 1:voxel.step_y * 2
            for i = 1:voxel.step_x * 2
                aa = aa + 1;
                mtrx2(j,i,k) = mtrx1{1,1}(aa,1);
            end
        end
    end
end

for k = 1:voxel.step_z * 2
    for j = 1:voxel.step_y * 2
        for i = 1:voxel.step_x * 2
            if mtrx2(j,i,k) == 1 || mtrx2(j,i,k) < max(mtrx1{1,1}) * SNR_filter % filter small SNR
                density_wt_cor(j,i,k) = 0;
                density(j,i,k) = 0;
            elseif rtio(j,i,k) < 0.1 % if there is a good enough SNR but the ratio is visibly small, don't count it
                    disp('!!! There in enough SNR in a suspiciosly fatty voxel !!!');
                    disp(strcat('ijk= ',num2str(i),',',num2str(j),',',num2str(k)));
                    density_wt_cor(j,i,k) = 0; %mtrx2(j,i,k);
                    density(j,i,k) =  0; %mtrx2(j,i,k) / rtio(j,i,k);
            else % associate the choline and density in each and compute the choline/density ratio:            
                    %density(j,i,k) = rtio(j,i,k);
                    density_wt_cor(j,i,k) = mtrx2(j,i,k);
                    density(j,i,k) = mtrx2(j,i,k) / rtio(j,i,k);
                    %density(j,i,k) = mtrx2(j,i,k);
            end            
        end
    end
end

% signal3D.wt_corr = interp3(X,Y,Z,density_wt_corr,XI,YI,ZI,'spline');
save(strcat(nfo1.PatientName.FamilyName,'_voxel.mat'),'voxel');
save(strcat(nfo1.PatientName.FamilyName,'_density.mat'),'density');
save(strcat(nfo1.PatientName.FamilyName,'_density_wt_cor'),'density_wt_cor');
save(strcat(nfo1.PatientName.FamilyName,'_rtio'),'rtio');


%% save the correlated data
ooo = 0;
for k = 1:c.pix_depth
    for j = 1:c.pix_width
        for i = 1:c.pix_height
            ooo = ooo + 1;
%            if abs(voxel.fov_cntr_x + round(shft_lr * voxel.FoV_x / voxel.notinterpfov_x) - ...
%        voxel.step_x * voxel.size_x) > abs(voxel.fov_cntr_x - voxel.p_fov_x / 2)
%                signal3D.list(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x) + 1;
%                signal3D.wt_cor(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x) + 1;
%            elseif abs(voxel.fov_cntr_x + round(shft_lr * voxel.FoV_x / voxel.notinterpfov_x) - ...
%        voxel.step_x * voxel.size_x) < abs(voxel.fov_cntr_x + voxel.p_fov_x / 2)
%                signal3D.list(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x) - 1;
%                signal3D.wt_cor(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x) - 1;
%            else
                signal3D.list(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x);
                signal3D.wt_cor(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x);
%            end
%            if abs(voxel.fov_cntr_y + round(shft_ud * voxel.FoV_y / voxel.notinterpfov_y) - ...
%        voxel.step_y * voxel.size_y) > abs(voxel.fov_cntr_y - voxel.p_fov_y / 2)
%                signal3D.list(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y) + 1;
%                signal3D.wt_cor(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y) + 1;
%            elseif abs(voxel.fov_cntr_y + round(shft_ud * voxel.FoV_y / voxel.notinterpfov_y) - ...
%        voxel.step_y * voxel.size_y) < abs(voxel.fov_cntr_y + voxel.p_fov_y / 2)
%                Ssignal3D.list(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y) - 1;
%                signal3D.wt_cor(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y) - 1;
%            else
                signal3D.list(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y);
                signal3D.wt_cor(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y);
%            end
            signal3D.list(ooo,3) = voxel.number_z / 2 + voxel.step_z - k + 1;
            signal3D.list(ooo,4) = density(j,i,k);
            signal3D.wt_cor(ooo,3) = voxel.number_z / 2 + voxel.step_z - k + 1;
            signal3D.wt_cor(ooo,4) = density_wt_cor(j,i,k);
        end
    end
end

signal3D.max = max(signal3D.list(:,4));
signal3D.mean = mean(signal3D.list(:,4));
dlmwrite('Corr_choSNR.txt', signal3D.list, 'delimiter', '\t', ...
         'precision', 6);
dlmwrite('Corr_choSNR_mean.txt', signal3D.mean, 'delimiter', '\t', ...
         'precision', 6);
dlmwrite('Corr_choSNR_max.txt', signal3D.max, 'delimiter', '\t', ...
         'precision', 6);
dlmwrite('Corr_choSNR_wt_corr.txt', signal3D.wt_cor, 'delimiter', '\t', ...
         'precision', 6);
%%
toc