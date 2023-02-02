function [lumen,epithelium_cells,epithelium,stroma] = gimmeSegs_prostate2(image,showMe,saveSegs,noDecon)
%% Prostate Adaptation of gimmeSegs

%% Purpose:  
% gimmeSegs modification for use with prostate slides
%
%% Inputs:
% image:  image being segmented
% showMe: enables figure of outputs
% saveSegs: saves outputs to .mat variables
% noDecon: use if not using Color Deconvolution
%
%% Outputs:
% Segmented lumen, epithelial cells, epithelium, stroma, and columnar cells

%% ------------------------------------------------------------------------
if ~exist('showMe','var')
    showMe = 0;
end

if ~exist('saveSegs','var')
    saveSegs = 0;
end

if ~exist('noDecon','var')
    noDecon = 0;
end

rrat = 0.293; % do not know what these are / ask Sam
brat = 0.5078;

disp('1. Segmenting Lumen')
tic

% lumen = imbinarize(image(:,:,2));
lumen = image(:,:,2)>190;
lumen = uint8(lumen);

toc
% -------------------------------------------------------------------------
disp('2. Segmenting Epithelium')

if noDecon == 1
    
    whatsLeft = image .* repmat(uint8(~lumen),[1,1,3]);
    whatsLeft = whatsLeft .* repmat(uint8(~vessel),[1,1,3]);
    
    red = whatsLeft(:,:,1);
    blue = whatsLeft(:,:,3);
    stroma = (red + blue) ./ 2;
    purp = ((red.*rrat) + (blue.*brat))./2;
    p2p = stroma - purp;
    
    nuc = uint8(zeros(size(p2p)));
    nuc(p2p>60) = 1;
    
    se = strel('disk',2);%nuc2 = imclose(nuc,se);
    nuc2 = imfill(nuc);
    cc = bwconncomp(nuc2);
    epithelium = uint8(zeros(size(nuc2)));
    x = regionprops('table',cc,'MajoraxisLength','MinoraxisLength','Circularity','Area');
    
    for i = 1:length(cc.PixelIdxList)
        if x.Circularity(i) > 0.1 && x.Area(i) > 1
            epithelium(cc.PixelIdxList{i}) = 1;
        end
    end
else
    % set of standard values for stain vectors (from python scikit)
    He = [0.6443186; 0.7166757; 0.26688856];
    Eo = [0.09283128; 0.9545457; 0.28324];
    Res = [ 0.63595444;   0.001; 0.7717266 ]; %residual
    
    HDABtoRGB = [He/norm(He) Eo/norm(Eo) Res/norm(Res)]';
    RGBtoHDAB = inv(HDABtoRGB);
    
    tic
    imageHDAB = SeparateStains(image, RGBtoHDAB);
    toc
    
    epithelium_cells = imcomplement(imageHDAB(:,:,1))>0.3; % >0.5
   
    thresh = quantile(unique(imcomplement(imageHDAB(:,:,1))),0.3);
    thresh2 = quantile(unique(imageHDAB(:,:,2)),0.6);
    epithelium = (imcomplement(imageHDAB(:,:,1))>(thresh+0.1)) + (imageHDAB(:,:,2)>thresh2) & ~lumen; %0.5 %0.8
    
%     epithelium2 = imbinarize(imcomplement(image(:,:,1)));
%     epithelium = (epithelium + epithelium2) & ~lumen;
    epithelium = uint8(epithelium);
    
%     figure('Position',[100 2000, 1500, 900]);
% subplot(221); imagesc(image); axis image
% subplot(222); imagesc(imoverlay(image,epithelium)); axis image

end

% -------------------------------------------------------------------------
disp('3. Segmenting Stroma')

tic

% stroma = imbinarize(imageHDAB(:,:,3));
thresh = quantile(unique(imcomplement(imageHDAB(:,:,2))),0.3); % find the 0.3 quantile to segment between epithelium and stroma
stroma = imcomplement(imageHDAB(:,:,2))>(thresh+0.1) & ~lumen; %0.4
% stroma2 = imageHDAB(:,:,1)>0.9; %0.8
% stroma = (stroma + stroma2) & ~lumen;
stroma = uint8(stroma);


%     figure('Position',[100 2000, 1500, 900]);
% subplot(221); imagesc(image); axis image
% subplot(222); imagesc(imoverlay(image,stroma)); axis image
toc

disp('Segmentation Complete');
if showMe == 1
    
    figure;
    subplot(141); imagesc(image); title('orig'); axis image
    subplot(142); imagesc(lumen); title('lumen'); axis image
    subplot(143); imagesc(epithelium); title('epithelium'); axis image
    subplot(144); imagesc(stroma); title('stroma'); axis image
    set(gcf,'Position',[100,100,1600,600]); axis image
    
end

if saveSegs == 1
    
    disp('Writing segmentations')
    
    save('lumen','lumen');
    save('epithelium','epithelium');
    save('stroma','stroma');
    
end