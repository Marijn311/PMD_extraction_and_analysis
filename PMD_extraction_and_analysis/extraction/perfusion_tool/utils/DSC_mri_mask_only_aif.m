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


function [mask]=DSC_mri_mask_only_aif(volumes,mask,options)


if options.display > 0
    disp(' ')
    disp('Using pre-calculated mask for the brain... ');
end



volume_sum=sum(volumes,4);

temp = mask;
clear mask;
mask.aif=temp;
hf_mask=zeros(options.nS,1);

for s=1:options.nS
    % fill any "holes" created by thresholding
    temp = mask.aif(:,:,s);
    
    % eliminate smaller connected components and keep those larger
    % than #options.mask.pixel
    
    CC = bwconncomp(temp);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    
    for j=1:size(numPixels,2)
        
        if numPixels(idx) > numPixels(j) && numPixels(j) < options.mask.npixel
            
            temp(CC.PixelIdxList{j}) = 0;
        end
    end
    
    mask.data(:,:,s)=imfill(temp,'holes');
    

end

mask.aif=mask.aif.*mask.data;



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

mask.aif = mask.data;
% Combine the initial mask with the refined mask
mask.aif=mask.aif.*mask.data;

end

