function msg = whole_mount_tiles_stitching
%% Whole Mount Tile Stitching

%% Purpose:  
% This function creates a whole-mount image from individual tiles and implements
% black_border_removal to remove border around image in the case of old
% images that were scanned with a border.

%% ------------------------------------------------------------------------

% directory = ''; name the needed directory that contains tiled images to
% be re-stitched together
% subj = XX; % patient number
% slice = "XX"; % slide name
% scanner = "Olympus"; % scanner name indictate tile sizes
% anot_micro = "Nikon"; % scanner used if originally annotated on a
% previous scanner

for d = 1:numel(directory)
    tic;
    addpath(directory(d))
    subj_slice = sprintf('%i_%s',subj(d),slice(d));
    fprintf('Processing: %s\n',subj_slice);
    
    disp('Loading images');
    
    files = fullfile(directory(d),'tile*.tiff');
    input_images = dir(files);
    input_images = {input_images.name}';
    
    if (isempty(input_images))
        files = fullfile(directory(d), 'tile*.jpg'); % adjust extension if needed
        input_images = dir(files);
        input_images = {input_images.name}';
    end
    
    name = strsplit(input_images{numel(input_images)},{'_','.'});
    xmax = str2num(name{2}(2:end));
    ymax = str2num(name{3}(2:end));
    
    x = [1 xmax];
    y = [1 ymax];
    
    xx = 1;
    yy = 1;
    
    
    if strcmp(scanner(d),'Huron')
        resize_factor = 10; % 8 nikon/olympus 12 huron
    else
        resize_factor = 8;
    end
    
    im = uint8(imread(input_images{numel(input_images)}));
    im_size = size(im);
    im = imrotate(imresize(im,1/8),90);
    ximage = cell(1,numel(input_images));
    yimage = cell(1,numel(input_images));
    for i = 1:numel(input_images)
        ximage{i} = ones(size(im(:,:,1)),'uint8') * xx;
        ximage{i} = imresize(ximage{i},[1 size(im,1)],'nearest');
        ximage{i} = imrotate(ximage{i},90);
        
        xx = xx + 1;
        if xx > xmax
            xx = 1;
        end
        
        yimage{i} = ones(size(im(:,:,1)),'uint8') * yy;
        yimage{i} = imresize(yimage{i},[1 size(im,2)],'nearest');
        yimage{i} = imrotate(yimage{i},90);
        
        yy = yy + 1;
        if yy > ymax
            yy = 1;
        end
    end
    
    disp('Creating rows');
    % Arrange rows
    ind1 = y(1);
    ind2 = y(2);
    rows = cell(1,xmax);
    for i = 1:xmax
        im = cell(1,ymax);
        ind3 = 1;
        for j = ind1:ind2
            im{ind3} = imrotate(imresize(uint8(imread(input_images{j})),1),90);
            ind3 = ind3 + 1;
        end

        rows{i} = imtile(im,'GridSize',[1 ymax]); clear im;
        rows{i} = imrotate(rows{i},180); % flips orientation
        
        ind1 = ind1 + ymax;
        ind2 = ind2 + ymax;
        if ind2 > numel(input_images)
            break
        end
    end
    
    toc
    xx = 1;
    xrows = cell(1,xmax);
    missing_vals = cell(1,xmax);
    for i = 1:xmax
        xrows{i} = ximage(xx:xmax:size(ximage,2)); % creates image rows
        xrows{i} = imtile(xrows{i},'GridSize',[1 size(rows{1},2)]);
        missing_vals{i} = uint8((xrows{i} == 0)*i); % xrows is filled with 0s beyond what was created above
        xrows{i} = xrows{i} + missing_vals{i};
        
        xx = xx + 1;
        if xx > xmax
            break
        end
    end
    
    yy = 1;
    ycols = cell(1,ymax);
    missing_vals = cell(1,ymax);
    for i = 1:ymax
        ycols{i} = yimage(yy:ymax:size(yimage,2)); % creates image cols
        ycols{i} = imtile(ycols{i},'GridSize',[1 xmax*size(rows{1},1)])'; % inverting creates col from row
        missing_vals{i} = uint8((ycols{i} == 0)*i);
        ycols{i} = ycols{i} + missing_vals{i};
        
        yy = yy + 1;
        
        if yy > ymax
            break
        end
    end
    
    toc
    
    disp('Compiling image');
       
    whole_image = imtile(rows,'GridSize',[xmax 1]);
    whole_image = imrotate(whole_image,-180);

    whole_ximage = flipud(imtile(xrows,'GridSize',[xmax 1])); % tile numbering starts in LL
    whole_ximage = double(whole_ximage);
    % figure(); imagesc(whole_ximage);
    
    whole_yimage = imtile(ycols,'GridSize',[1 ymax]);
    whole_yimage = double(whole_yimage);
    % figure(); imagesc(whole_yimage);
    toc
    disp('Saving image');
%     figure(); imshow(whole_image);
    %                 saveas(gcf,sprintf('%s_whole_mount_recon_%i',subj_slice,resize_factor),'tif');

toc
msg = 'done';
disp(msg)
