% minarikova.lenka@gmail.com

function [signal3D] = psf_pics(directory, no_cor, fctr, CSI_shft_ud, CSI_shft_lr)

% directory = '/Users/zgung/Desktop/HA/';
% no_cor = 0;
% fctr = 2;
% CSI_shft_ud = 0;
% CSI_shft_lr = 0;

% !!!!!!!!!!!!!!!!!! readme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% for working you need my other function called read_ascconv_lenk.m
%   for reading parameters from the dicom
%   and you need a txt file with a set of SNR to evaluate, made with
%   my other function called SNR_lenk.m
% directory = '~/Patient_name' - where is a directory called "Spec" with
%   a dicom 3D CSI file
% fctr = 2 for patients, 4 for phantoms with 2 times higher resolution in CSI
% no_cor = 0 for corrected maps, 1 for initial (not corrected) maps, 2 for ratio



% the output is maximal, mean value of all SNRs of Cho and a table
%   with all SNRs in one row, all saved in txt files in Spec directory

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% search for variables in spectroscopy file:
tic;

%directory_main = directory;
directory_spec = strcat(directory,'Spec/');
cd(directory_spec);
disp(strcat('Processing:',directory));
%% <<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>
%csishift = 0; %if the matrix is shifted of 25% (4.2) up

% read from spectroscopic file from the text part
slices = dir(pwd);
for k = length(slices):-1:1 % find .IMA file
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    elseif strcat(fname(end-2),fname(end-1),fname(end)) ~= strcat('IMA')
        slices(k) = [];
    end
end
%% check if the slices are in a good order
Afields = fieldnames(slices);
Acell = struct2cell(slices);
sz = size(Acell);
% Convert to a matrix
Acell = reshape(Acell, sz(1), []);      % Px(MxN)
% Make each field a column
Acell = Acell';                       % (MxN)xP
k = 1;
for k = 1:length(slices)
    [~,~,~,~,slc_n,~,~,~,~,~,~,~,~,~] = strread(slices(k,1).name,'%s %s %s %d %d %d %d %d %d %d %d %d %d %s','delimiter','.');
    %slc_n = num2cell(slc_n);
    Acell{k,4}= slc_n;
end
% Sort by first field "name"
Acell = sortrows(Acell, 4);
% Put back into original cell array format
Acell = reshape(Acell', sz);
% Convert to Struct
slices = cell2struct(Acell, Afields, 1);
%voxel = read_ascconv_lenk(slices(1).name); % read parrameter from spectroscopy dicom
nfo1 = dicominfo(slices(1).name);
%[signal3D,sgmnts,voxel] = psf(slices(1).name);
%%

%load(strcat(nfo1.PatientName.FamilyName,'_sgmnts_w.mat')); % load the segmented water data
load(strcat(nfo1.PatientName.FamilyName,'_voxel.mat'));
load(strcat(nfo1.PatientName.FamilyName,'_density.mat'));
load(strcat(nfo1.PatientName.FamilyName,'_imgs_w.mat'));
load(strcat(nfo1.PatientName.FamilyName,'_density_wt_cor'));
load(strcat(nfo1.PatientName.FamilyName,'_rtio'));
%%
c.cut_width = round(2 * voxel.step_y * voxel.size_y); % width in pixels of the cuted image
c.cut_height = round(2 * voxel.step_x * voxel.size_x); % height -||-
c.cut_depth = round(2 * voxel.step_z * voxel.size_z); % depth -||-
c.pix_width = voxel.step_y * 2; % number of voxels in the pressbox - x axis
c.pix_height = voxel.step_x * 2; % -||- - y axis
c.pix_depth = voxel.step_z * 2; % -||- - z axis

if no_cor == 1
    density = density_wt_cor;
elseif no_cor == 2
    density = rtio;
end
%%
f_nas = fctr * voxel.size_x;
c.density = size(density);
[XI,YI,ZI] = meshgrid(1:(c.density(2) - 0.99) / (c.density(2) * f_nas):c.density(2),...
    1:(c.density(1) - 0.99) / (c.density(1) * f_nas):c.density(1),...
    1:(c.density(3) - 0.99) / (c.density(3) * f_nas):c.density(3));
%%
signal3D_ = ba_interp3(density,XI,YI,ZI,'linear');
clear XI; clear YI; clear ZI;
%% interpolate images so they have "fctr" times bigger resolution
c.img_s = size(imgs_w);
[XI,YI,ZI] = meshgrid(1:(c.img_s(2) - 0.9) / (c.img_s(2) * fctr):c.img_s(2),...
    1:(c.img_s(1) - 0.9) / (c.img_s(1) * fctr):c.img_s(1),...
    1:(c.img_s(3) - 0.9) / (c.img_s(3) * fctr):c.img_s(3));
imgs_w = ba_interp3(imgs_w,XI,YI,ZI,'linear');
clear XI; clear YI; clear ZI;
%% define the colormaps for background MRI picture and metabolit map
% than find a maximal value of background and signal images:
h = max(max(max(imgs_w)));
g_test = max(max(max(signal3D_)));
%%
if g_test < 10
    g = max(max(max(signal3D_))) * 100; % amplified signal values
else
    g = g_test;
end
%mx = max(h);
cmp = colormap(gray);
[X,Y] = meshgrid(1:3,1:64);
[XI,YI] = meshgrid(1:3,1:64 / h:64);
cmp1 = interp2(X,Y,cmp,XI,YI);
jet1 = colormap(jet);
[X,Y] = meshgrid(1:3,1:64);
[XI,YI] = meshgrid(1:3,1:64 / g:64);
jet2 = interp2(X,Y,jet1,XI,YI);

%% draw 8 slices at once:
figure('units','normalized','position',[.05 .1 .95 .95])
for aa = 1:c.pix_depth
    %subplot('Position',[(-0.0 ) 0.05 (2.5 / c.pix_depth) (2.5 / c.pix_depth)]);
    subaxis(round(c.pix_depth / 3) + 1,5 ,aa, 'Spacing', 0, 'Padding', 0, 'Margin', 0);
    subimage(imgs_w(:,:,fix((aa) * (c.cut_depth * fctr) / (c.pix_depth + 1)) + 5),cmp1);
    pic = signal3D_(:,:,fix((aa) * (c.cut_depth * fctr) / (c.pix_depth + 1)) - 0);
    hold on
    pic_norm = pic/max(max(pic));
%     if g_test < 10
%         hh2 = image((voxel.fov_x1 - 3) * fctr,(voxel.fov_y1 - 3) * fctr,pic .* 100,...
%             'AlphaDataMapping','scaled',...
%             'AlphaData',gradient(pic));
%     else
        hh2 = image((voxel.fov_x1 - 3) * fctr,(voxel.fov_y1 - 3) * fctr,pic); %,...
%            'AlphaDataMapping','scaled',...
%            'AlphaData',pic/15);
%     end
    
    set( hh2, 'AlphaData', .3);
    colormap(jet2);
    axis([100 800 100 800]);
    %axis([0 length(imgs_w{1,1}) 0 numel(imgs_w{1,1}(:,1))]);
    caxis([(0 * g) (1 * g)]);
    hold off
    axis off
end

toc