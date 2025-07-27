% Calculates masks for DSC-MRI exams.
%
% Inputs: 
% - volumes: (4D Matrix) containing the DSC signal of all voxels over time.
% - options: struct containing the method options, the significant ones are:
%   - options.mask.npixel: represents the minimum number of pixels of a connected component
%                          used as a threshold to exclude the scalp and areas adjacent to the outside of the brain.
%   - options.display: 
%     - level 1: Shows the progress of the processing.
%     - level 2: Shows the masks and information on the threshold and on the intensities of the images to be masked.
%
% Outputs: 
% - mask: structure containing:
%   - aif: optimized mask for the search of the arterial input function.
%   - data: optimized mask for masking the entire brain.
%   - figure 2: an example image sliceS showing the masked brain regions. 
%               The mask edges are highlighted in green, and the original intensity data is shown in grayscale. 
%               This visualization helps in verifying the accuracy of the mask and ensuring that 
%               the brain regions are correctly identified and segmented from the surrounding tissues.


function [mask]=DSC_mri_mask(volumes,options)

% Display progress if options.display > 0
if options.display > 0
    disp(' ')
    disp('Masking data... ');
end
 
% Sum the volumes across the 4th dimension (time) to get a single 3D volume
volume_sum=sum(volumes,4);
mask.data=zeros(size(volume_sum)); % Initialize mask for brain data
mask.aif=zeros(size(volume_sum)); % Initialize mask for arterial input function (AIF)

% Make a 2D binary mask, this is needed as input for the active contours. 
% Extracts the 10th slice (2D image) from the 3D volume DSC_volume and creates a binary mask (seed). 
% The mask is created by setting elements to true (1) where the corresponding elements in the slice are greater than 75, and false (0) otherwise.

%initialize the 3D seed mask
seed_mask = volume_sum() > 100;

% Perform the segmentation using active contours, specifying the seed mask.
mask.data = activecontour(volume_sum, seed_mask,300);
mask.data = double(mask.data); % Convert mask from boolean to numeric (0 or 1)

% Per slice: Fill the holes in the mask, and remove small unconnected components 
for s=1:options.nS
    mask_slice = mask.data(:,:,s); % take a slice of the mask
    
    % Sometimes a bit of the skull is still attach in the segmentation because some part of the skull was very close to the brain.
    % To remove this, we erode the mask a bit. This will hopefully un-connect the skull fragments from the brain in 3D.
    % Then we remove the all components that are not attached to the main component (which is the brain).  
    % Then we dilate the mask again to get the original size of the brain mask back.
    % This also smooths the mask a bit.
    structure_element = strel('disk', 5);
    mask_slice = imerode(mask_slice, structure_element);

    % Eliminate smaller connected components and keep those larger than options.mask.npixel
    CC = bwconncomp(mask_slice);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    
    for j=1:size(numPixels,2)
        if numPixels(idx) > numPixels(j) && numPixels(j) < options.mask.npixel
            mask_slice(CC.PixelIdxList{j}) = 0;
        end
    end
    
    mask_slice = imdilate(mask_slice, structure_element);

    % Fill holes in the slice and update the mask to new include the slie version with filled holes
    mask.data(:,:,s)=imfill(mask_slice,'holes');
    
    
end




% Remove all components that are not attached to the main component in 3d
CC = bwconncomp(mask.data);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
for j=1:size(numPixels,2)
    if numPixels(idx) > numPixels(j)
        mask.data(CC.PixelIdxList{j}) = 0;
    end
end






if options.display > 0
    disp('Active contour mask generated successfully');
end


% Create a figure with a slider to scroll through slices
if options.display > 1
    middle_slice_index = round(size(volume_sum,3)/2);
    hf_slider = figure('Name', 'Slice Viewer');
    slice_slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', options.nS, 'Value', middle_slice_index, ...
                                'Units', 'normalized', 'Position', [0.1 0.01 0.8 0.05], ...
                                'Callback', @update_slice);
    slice_text = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.45 0.07 0.1 0.05], ...
                            'String', ['Slice: ' num2str(middle_slice_index)]);
    ax = axes('Position', [0.1 0.2 0.8 0.7]);
    update_slice(slice_slider);
end


function update_slice(src, ~)
    slice_idx = round(get(src, 'Value'));
    set(slice_text, 'String', ['Slice: ' num2str(slice_idx)]);
    imagesc(ax, volume_sum(:,:,slice_idx));
    colormap(ax, 'gray');
    axis(ax, 'image');
    title(ax, ['Slice ' num2str(slice_idx) '/' num2str(options.nS) ' - Example of Masked Data']);
    impixelinfo;


    B = bwboundaries(mask.data(:,:,slice_idx), 4);
    hold(ax, 'on');
    for k = 1:length(B)
        boundary = B{k};
        plot(ax, boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
    end
    hold(ax, 'off');
end


%i will change the aif selection so for now the aif mask is the same as the data mask
mask.aif = mask.data;
% Combine the initial mask with the refined mask
mask.aif=mask.aif.*mask.data;

% Close figures only when the user presses the 'x' in the popup window
if options.display > 0
    if options.display > 2
        for cont_h = 1:length(hf_mask)
            set(hf_mask(cont_h), 'CloseRequestFcn', @my_closefcn);
        end
    end
end

function my_closefcn(src, ~)
    delete(src);
end

end

