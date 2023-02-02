function image_no_border = black_border_removal(image)
%% Black Border Removal

%% Purpose:  
% This function finds the maximum and minimum image rows/cols to remove
% border around whole slides images.
%
%% Inputs:
% image: image with a border
%
%% Outputs:
% image_no_border: image without border

%% ------------------------------------------------------------------------
% Remove black border
% if ref > 0 
   grayImage = rgb2gray(image);
   grayImage(find(grayImage==255)) = 0; % replace white border with black
    [nonZeroRows,nonZeroColumns] = find(grayImage); % find nonzero rows/cols

    % Get the cropping parameters
    topRow = min(nonZeroRows(:));
    bottomRow = max(nonZeroRows(:));
    leftColumn = max(nonZeroColumns(:));
    rightColumn = min(nonZeroColumns(:));

    crop_parameters = [rightColumn topRow leftColumn bottomRow];
    image_no_border = imcrop(image,crop_parameters);
% else
%     grayImage = im2gray(ref);
%     [nonZeroRows,nonZeroColumns] = find(grayImage); % find nonzero rows/cols
% 
%     % Get the cropping parameters
%     topRow = min(nonZeroRows(:));
%     bottomRow = max(nonZeroRows(:));
%     leftColumn = max(nonZeroColumns(:));
%     rightColumn = min(nonZeroColumns(:));
% 
%     crop_parameters = [rightColumn topRow leftColumn bottomRow];
%     image_no_border = imcrop(image,crop_parameters);

% end 
end