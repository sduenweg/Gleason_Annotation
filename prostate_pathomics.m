 function status = prostate_pathomics
% Pathomic Feature Calculator

%% ------------------------------------------------------------------------
show_results = 0;

fid = fopen(sprintf('tile_seg_output1-%s.csv',datestr(now(),'mm_dd_yyyy')),'a');
fid2 = fopen(sprintf('tile_seg_output2-%s.csv',datestr(now(),'mm_dd_yyyy')),'a');

anot = {'Seminal_Vesicles','Atrophy','HGPIN','G3','G4FG','G4CG','G5','Tissue'}; % annotation classes
scanner = {'Huron','Olympus'}; % whole slide image scanners

for s = 1:numel(scanner)
    d = dir(sprintf('%s/*/*.tiff',scanner{s}));
    date = {d.date}';

    input_images = {d.name}';
    folder = {d.folder}';
    
    for i = 1:numel(d)
        grade_dir = strsplit(folder{i},'/');
        grade{i} = grade_dir{end};
        grade_number(i) = find(contains(anot,grade{i}));
    end
    
    for i = 1:numel(d)
        grade_dir = strsplit(folder{i},'/');
        grade{i} = grade_dir{end};
        grade_number(i) = find(contains(anot,grade{i}));
    end
    
    for x = 1:numel(input_images)
        fprintf('Processing: %s\n',input_images{x});
        im = imread(input_images{x});
        im_resized = imresize(im,0.4);
        
        % Image color deconvolution
        im_resized = im_resized(:,:,1:3);
        I = rgb2gray(im_resized);
        
        tissue_mask = I<200;
        tissue_mask = ~bwareaopen(~tissue_mask,6000); % removes rip noise
        tissue_mask = imgaussfilt(uint8(tissue_mask),4);
        tissue_mask = uint8(bwareaopen(tissue_mask,50000)); % help flood lumen
        tissue_mask_label = bwlabel(tissue_mask);
        
        [lumen_mask,epithelium_cells,epithelium_mask,stroma_mask] = gimmeSegs_prostate2(im_resized,0,0,0);
        
        disp('Cleaning up segmentations');
        epithelium_mask = epithelium_mask .* tissue_mask;
        stroma_mask = stroma_mask .* tissue_mask;
        lumen_mask_noise = lumen_mask .* tissue_mask; % masks small lumen and some epithelium
        lumen_mask_noise = bwareafilt(logical(lumen_mask_noise),[0 150]);
        
        % Clean up lumen
        if strcmp(scanner{s},'Huron')
            lumen_mask = im_resized(:,:,2)>240; % huron lumen more white than oly
        end
        se = strel('disk',1);
        lumen_mask = lumen_mask & ~lumen_mask_noise;
        lumen_mask = imclose(lumen_mask,se);
        lumen_mask = imfill(lumen_mask,'holes');
        lumen_mask = bwareaopen(lumen_mask,100); %100
        lumen_mask = imgaussfilt(uint8(lumen_mask),1);
% -------------------------------------------------------------------------
        % Clean up stroma
        stroma_mask = stroma_mask & ~lumen_mask & ~epithelium_cells;
        stroma_mask = bwmorph(stroma_mask,'bridge');
        stroma_mask = ~bwareaopen(~stroma_mask,200);
        stroma_mask = bwareaopen(stroma_mask,20);
        stroma_mask = imgaussfilt(uint8(stroma_mask),1);
        
% -------------------------------------------------------------------------
        % Clean up epithelium
        epithelium_mask = logical(epithelium_mask) + epithelium_cells;
        
        if grade_number(x) < 4 | grade_number(x) == 8
            se = strel('disk',4);
            epithelium_mask = imclose(epithelium_mask,se);
        elseif grade_number(x) == 6
            se = strel('disk',1);
            epithelium_mask = imclose(epithelium_mask,se);
        end
        
        epithelium_mask = logical(epithelium_mask) & ~lumen_mask & ~stroma_mask;
        
        epithelium_mask = bwmorph(epithelium_mask,'thick',2);
        epithelium_mask = bwareaopen(epithelium_mask,250);
        epithelium_mask = ~bwareaopen(~epithelium_mask,500);
        epithelium_mask = bwareaopen(epithelium_mask,500);
        epithelium_mask = imgaussfilt(uint8(epithelium_mask),1);
        
        epithelium_mask_filled = imfill(epithelium_mask,'holes') & ~lumen_mask & ~stroma_mask;
        epithelium_mask_filled = ~bwareaopen(~epithelium_mask_filled,100);
        epithelium_mask_filled = bwareaopen(epithelium_mask_filled,500);
        
        epithelium_mask_filled = imgaussfilt(uint8(epithelium_mask_filled),1);
        
% -------------------------------------------------------------------------
        lumen_mask = bwmorph(lumen_mask,'thick',2);
        lumen_mask = uint8(imclearborder(lumen_mask,4)); % removes any lumen touching tile edges
        
        disp('Segmenting Gland');
        SEL_image = stroma_mask + 2*lumen_mask + 3*epithelium_mask_filled;
        
        %         figure('Position',[100 2000, 1500, 900]);
        %         subplot(131); imagesc(im_resized); axis image
        %         title(input_images{x});
        %         subplot(132); imagesc(SEL_image); axis image
        %         title('SEL_image');
        %         subplot(133); imagesc(SEL_image_epi); axis image
        %         title('SEL_image separate epi');
        %         saveas(gcf,sprintf('SEL_%s',input_images{x}),'tif');
        
        
        % Begin epith_wall_thickness_calculation, so re-map the variables into the legacy code
        bw_lumen = uint8(lumen_mask);
        bw_lumen_label = bwlabel(bw_lumen);
        bw_epith_label = bwlabel(epithelium_mask);
        bw_stroma = stroma_mask;
        bw_epith = uint8(epithelium_mask);
        bw_epith_filled = epithelium_mask_filled;
        bw_epith_filled_label = bwlabel(bw_epith_filled);
        core_mask = tissue_mask;
        bw_cells = uint8(epithelium_cells);
        
        %         figure('Position',[100 2000, 1500, 900]);
        %         subplot(2,3,1);
        %         imagesc(im_resized); axis image;
        %         title('Raw Image')
        %         subplot(2,3,2);
        %         imagesc(tissue_mask_label); axis image;
        %         title('Tissue Mask')
        %         subplot(2,3,3);
        %         imagesc(bw_epith);axis image;
        %         title('Epithelium Mask');
        %         subplot(2,3,4);
        %         imagesc(bw_epith_label);axis image;
        %         title('Epithelium labels');
        %         subplot(2,3,5);
        %         imagesc(bw_lumen);axis image;
        %         title('Lumen Mask')
        %         subplot(2,3,6);
        %         imagesc(bw_stroma);axis image;
        %         title('Stroma Mask')
        %         % saveas(gcf,sprintf('seg_%s',input_images{x}),'tif');
        
% -------------------------------------------------------------------------
        % Stroma v epithelium area calculations
        stroma_area = bwarea(bw_stroma);
        epith_area = bwarea(bw_epith);
% -------------------------------------------------------------------------
        
        [Bl,Ll,Nl,Al] = bwboundaries(bw_lumen_label,'noholes');
        statsl = regionprops(Ll,'Area','Centroid','PixelList');
        
        [Be,Le,Ne,Ae] = bwboundaries(bw_epith_filled_label,'noholes');
        statse = regionprops(Le,'Area','Centroid','PixelList');
        
        colors=['b' 'g' 'r' 'c' 'm' 'y'];
        %         colors = ['c' 'c' 'c' 'c' 'c' 'c'];
        
        ep_lb = [];
        
%% ------------------------------------------------------------------------
        bw_lumen_roundness = bw_lumen_label;
        bw_lumen_area = bw_lumen_label;
        bw_epith_thickness = bw_epith_filled_label;
        bw_cell_frac = bw_epith_filled_label;
        bw_epith_size = bw_epith_filled_label;
        bw_epith_roundness = bw_epith_filled_label;
        lumen_tort = [];
        epith_size = [];
        cell_frac = [];
        wall_thickness = [];
        min_lum_index_thick = [];
        wall_thick = [];
        epith_tort = [];
        
%         figure('Position',[100 100 1500 1500])
%         imagesc(im_resized);axis image; hold on;
        
        for i = 1:length(Bl)
            ep_lb(i) = bw_epith_filled_label(statsl(i).PixelList(1,2),statsl(i).PixelList(1,1)); % this maybe was backwards? 2,1?
            if ep_lb(i) == 0 %This indicates that an epithelium didn't get closed, and these are then excluded
                disp(sprintf('lumen number %i is connected with an epithelium edge and will be excluded',i));
                
                j = i; % in next step j is used for loop
                boundaryl = Bl{j};
                area = statsl(j).Area;
                delta_sq = diff(boundaryl).^2;
                perimeter = sum(sqrt(sum(delta_sq,2)));
                lumen_tort(j) = 4*pi*area/(perimeter^2);
                
                bw_lumen_roundness(bw_lumen_roundness == j) = lumen_tort(j); % replaces labeled lumen with tortuosity value
                bw_lumen_area(bw_lumen_area == j) = area;
                
                try
                    ep_lb(i) = bw_epith_filled_label(statsl(i).PixelList(2,2),statsl(i).PixelList(2,1));
                catch
                    disp('hey-you');
                end
            end
        end
        
        j = 0;
        inde = 0;
        for i = 1:length(Be) % Loop over the epithelium
            num_skipped = 0;
            if statse(i).Area<100000000 && statse(i).Area>0
                % plot the epithelium
                inde = inde+1;
                output_data.data(x).gland(inde).boundarye = Be{i};
                %             output_data.data(x).gland(inde).epith_size = statse(i).Area;
                ep_masked_cells_im = uint8(bw_epith_filled_label == i).* bw_cells;
                ep_masked_lumen_im = uint8(bw_epith_filled_label == i).* bw_lumen;
                output_data.data(x).gland(inde).epith_size = statse(i).Area - numel(find(ep_masked_lumen_im));
                output_data.data(x).gland(inde).cell_frac = sum(sum(ep_masked_cells_im)) / (statse(i).Area - sum(sum(ep_masked_lumen_im)));
                boundarye = Be{i};
                ep_com = [round(mean(boundarye(:,1))) round(mean(boundarye(:,2)))];
                %             current_core = tissue_mask_label(ep_com(1),ep_com(2));
                skip_pic = 0;
                try
                    glnd_im_out = im_resized(ep_com(1)-111:ep_com(1)+112,ep_com(2)-111:ep_com(2)+112,:);
                catch
                    skip_pic = 1;
                end
                cidx = mod(i,length(colors))+1;
                plot(boundarye(:,2), boundarye(:,1),...
                    colors(cidx),'LineWidth',2);
                
                areae = statse(i).Area;
                delta_sqe = diff(boundarye).^2;
                perimetere = sum(sqrt(sum(delta_sqe,2)));
                output_data.data(x).gland(inde).epith_tort = 4*pi*areae/perimetere^2;
                
                % For plots -----------------------------------------------
                epith_size(i) = statse(i).Area - numel(find(ep_masked_lumen_im));
                cell_frac(i) = sum(sum(ep_masked_cells_im)) / (statse(i).Area - sum(sum(ep_masked_lumen_im)));
                epith_tort(i) = 4*pi*areae/perimetere^2;
                %----------------------------------------------------------
                
                %randomize text position for better visibility
                rndRow = ceil(length(boundarye)/(mod(rand*i,7)+1));
                col = boundarye(rndRow,2); row = boundarye(rndRow,1);
                
                % Figure out the lumen within
                indl = 0;
                if numel(find(ep_lb==i))>0
                    % loop over the lumen within:
                    for j = find(ep_lb==i)
                        try
                            
                            if statsl(j).Area>1
                                indl = indl + 1;
                                % plot the lumen boundaries
                                boundaryl = Bl{j};
                                cidx = mod(i,length(colors))+1;
                                plot(boundaryl(:,2), boundaryl(:,1),...
                                    colors(cidx),'LineWidth',2);
                                %randomize text position for better visibility
                                rndRow = ceil(length(boundaryl)/(mod(rand*i,7)+1));
                                col = boundaryl(rndRow,2); row = boundaryl(rndRow,1);
                                
                                area = statsl(j).Area;
                                delta_sq = diff(boundaryl).^2;
                                perimeter = sum(sqrt(sum(delta_sq,2)));
                                output_data.data(x).gland(inde).lumen(indl).lumen_tort = 4*pi*area/perimeter^2;
                                output_data.data(x).gland(inde).lumen(indl).wall_thick = 0;
                                output_data.data(x).gland(inde).lumen(indl).area = area;
                                
                                % For plots -------------------------------
                                lumen_tort(j) = 4*pi*area/perimeter^2;
                                bw_lumen_roundness(bw_lumen_roundness == j) = lumen_tort(j); % replaces labeled lumen with tortuosity value
                                bw_lumen_area(bw_lumen_area == j) = area;
                                % -----------------------------------------
                                
                                wall_thickness = [];
                                min_lum_index_thick = [];
                                for jj = 1:size(Bl{j},1)
                                    ind = 1;
                                    for k = 1:5:size(Be{i},1)
                                        wall_thickness(ind) = sqrt((Bl{j}(jj,1)-Be{i}(k,1))^2 + (Bl{j}(jj,2)-Be{i}(k,2))^2);
                                        ind = ind+1;
                                    end
                                    min_lum_index_thick(jj) = min(wall_thickness);
                                    if show_results == 1
                                        if rem(jj,30)==0
                                            h = text(Bl{j}(jj,2), Bl{j}(jj,1), sprintf('%0.2f',min_lum_index_thick(jj)));
                                            set(h,'Color',colors(cidx),...
                                                'FontSize',8);
                                        end
                                    end
                                end
                                
                                output_data.data(x).gland(inde).lumen(indl).wall_thick = mean(min_lum_index_thick);
                                wall_thick(i) = mean(min_lum_index_thick);
                                
                                if show_results == 1
                                    h = text(Bl{j}(1,2)+20, Bl{j}(1,1)+10, sprintf('%i\n%i\n%0.2f\n%0.2f',statse(i).Area,statsl(j).Area,output_data.data(x).gland(inde).lumen(indl).lumen_tort,output_data.data(x).gland(inde).lumen(indl).wall_thick))  ;
                                    set(h,'Color',colors(cidx),...
                                        'FontSize',8);
                                end
                                
                                % output_data.data(x).gland(inde)
                                % disp('results shown');
                                %
%                                 try
                                fprintf(fid,'%s,%s,%i,%s,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f\n',...
                                    scanner{s},input_images{x},x,grade{x},...
                                    grade_number(x),inde,indl,stroma_area,epith_area,...
                                    output_data.data(x).gland(inde).epith_size,...
                                    output_data.data(x).gland(inde).epith_tort,...
                                    output_data.data(x).gland(inde).cell_frac,...
                                    output_data.data(x).gland(inde).lumen(indl).area,...
                                    output_data.data(x).gland(inde).lumen(indl).wall_thick,...
                                    output_data.data(x).gland(inde).lumen(indl).lumen_tort);
                            end
                        catch
                            disp('Debug point')
                        end
                    end
                    
                    try
                        
                        tot_lum_area = [];
                        tot_wall_thick = [];
                        tot_lum_tort = [];
                        
                        for mmm = 1:numel(output_data.data(x).gland(inde).lumen(:))
                            tot_lum_area = [tot_lum_area output_data.data(x).gland(inde).lumen(mmm).area];
                            tot_wall_thick = [tot_wall_thick output_data.data(x).gland(inde).lumen(mmm).wall_thick];
                            tot_lum_tort = [tot_lum_tort output_data.data(x).gland(inde).lumen(mmm).lumen_tort];
                        end
                        fprintf(fid2,'%s,%s,%i,%s,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f\n',...
                            scanner{s},input_images{x},x,grade{x},grade_number(x),inde,...
                            output_data.data(x).gland(inde).epith_size,...
                            output_data.data(x).gland(inde).epith_tort,...
                            output_data.data(x).gland(inde).cell_frac,...
                            mean(tot_lum_area),...
                            mean(tot_wall_thick),...
                            mean(tot_lum_tort),...
                            sum(sum(bw_lumen))/sum(sum(core_mask)),...
                            sum(sum(bw_epith))/sum(sum(core_mask)),...
                            sum(sum(bw_stroma))/sum(sum(core_mask)));
                        
                    catch
                        disp('debug point');
                    end
                    
                else
                    disp('skipping epith, no lumen');
                    
                    ep_masked_cells_im = uint8(bw_epith_filled_label == i).* bw_cells;
                    ep_masked_lumen_im = uint8(bw_epith_filled_label == i).* bw_lumen;
                    epith_size(i) = statse(i).Area - numel(find(ep_masked_lumen_im));
                    cell_frac(i) = sum(sum(ep_masked_cells_im)) / (statse(i).Area - sum(sum(ep_masked_lumen_im)));
                    
                    boundarye = Be{i};
                    areae = statse(i).Area;
                    delta_sqe = diff(boundarye).^2;
                    perimetere = sum(sqrt(sum(delta_sqe,2)));
                    epith_tort(i) = 4*pi*areae/perimetere^2;
                    
                end
                
            else
                disp('skipping epith, too big');
                
                boundaryl = Bl{j};
                area = statsl(j).Area;
                delta_sq = diff(boundaryl).^2;
                perimeter = sum(sqrt(sum(delta_sq,2)));
                lumen_tort(j) = 4*pi*area/(perimeter^2);
                
                bw_lumen_roundness(bw_lumen_roundness == j) = lumen_tort(j); % replaces labeled lumen with tortuosity value
                bw_lumen_area(bw_lumen_area == j) = area;
                
                ep_masked_cells_im = uint8(bw_epith_filled_label == i).* bw_cells;
                ep_masked_lumen_im = uint8(bw_epith_filled_label == i).* bw_lumen;
                epith_size(i) = statse(i).Area - numel(find(ep_masked_lumen_im));
                cell_frac(i) = sum(sum(ep_masked_cells_im)) / (statse(i).Area - sum(sum(ep_masked_lumen_im)));
                
                boundarye = Be{i};
                areae = statse(i).Area;
                delta_sqe = diff(boundarye).^2;
                perimetere = sum(sqrt(sum(delta_sqe,2)));
                epith_tort(i) = 4*pi*areae/perimetere^2;
                
            end
%% -------------------------------------------------------------
            % display epithelium tort, size, and cell frac
            bw_epith_roundness(bw_epith_roundness == i) = epith_tort(i);
            bw_epith_size(bw_epith_size == i) = epith_size(i);
            bw_cell_frac(bw_cell_frac == i) = cell_frac(i);
            try
                bw_epith_thickness(bw_epith_thickness == i) = wall_thick(i);
            catch
                bw_epith_thickness(bw_epith_thickness == i) = 0;
            end
%% -------------------------------------------------------------
        end
%% ------------------------------------------------------------------------
%         % display features
%         figure('Position',[100 100 1500 1500])
%         imshow(im_resized); axis image; hold on;
%         h = imagesc(bw_lumen_roundness); caxis([0 0.7]); % change caxis to see colormap better
%         hold off
%         
%         [M,N] = size(bw_lumen_roundness);
%         alpha_data = bw_lumen_roundness > 0;
%         alpha_data = alpha_data(1:M, 1:N);
%         set(h, 'AlphaData', alpha_data);
%         title('lumen tort 2');
%         
%         figure('Position',[100 100 1500 1500])
%         imshow(im_resized); axis image; hold on;
%         h = imagesc(bw_lumen_area); caxis([0 3*10^4]); % change caxis to see colormap better
%         hold off
%         
%         [M,N] = size(bw_lumen_area);
%         alpha_data = bw_lumen_area > 0;
%         alpha_data = alpha_data(1:M, 1:N);
%         set(h, 'AlphaData', alpha_data);
%         title('lumen area 2');
%         
%         figure('Position',[100 100 1500 1500])
%         imshow(im_resized); axis image; hold on;
%         h = imagesc(bw_epith_thickness);
%         hold off
%         
%         [M,N] = size(bw_epith_thickness);
%         alpha_data = bw_epith_thickness > 0;
%         alpha_data = alpha_data(1:M, 1:N);
%         set(h, 'AlphaData', alpha_data);
%         title('epith thick 2')
%         
%         figure('Position',[100 100 1500 1500])
%         imshow(im_resized); axis image; hold on;
%         h = imagesc(bw_epith_roundness); caxis([0 1]); % change caxis to see colormap better
%         hold off
%         
%         [M,N] = size(bw_epith_roundness);
%         alpha_data = bw_epith_roundness > 0;
%         alpha_data = alpha_data(1:M, 1:N);
%         set(h, 'AlphaData', alpha_data);
%         title('epith tort 2')
%         
%         figure('Position',[100 100 1500 1500])
%         imshow(im_resized); axis image; hold on;
%         h = imagesc(bw_epith_size); %caxis([0 250]);
%         hold off
%         
%         [M,N] = size(bw_epith_size);
%         alpha_data = bw_epith_size > 0;
%         alpha_data = alpha_data(1:M, 1:N);
%         set(h, 'AlphaData', alpha_data);
%         title('epith_size 2');
%         
%         figure('Position',[100 100 1500 1500])
%         imshow(im_resized); axis image; hold on;
%         h = imagesc(bw_cell_frac); caxis([0 1]);
%         hold off
%         
%         [M,N] = size(bw_cell_frac);
%         alpha_data = bw_cell_frac > 0;
%         alpha_data = alpha_data(1:M, 1:N);
%         set(h, 'AlphaData', alpha_data);
%         title('cell frac 2')
%% ------------------------------------------------------------------------
        fprintf('%s done\n',input_images{x});
        % saveas(gcf,sprintf('outline_seg_%s',input_images{x}),'tif');
        close all
        
    end
end
fclose(fid);
status = 'done';
end