% Function of the DSC_mri package - DSC_mri_conc
%
% Calculates the contrast agent concentrations and S0 maps in DSC-MRI exams.
%
% Input parameters:
% - volumes (4D matrix) containing the DSC signal trends of all voxels.
% - mask (3D matrix) contains the matrix to mask every not belonging to the brain   interest for the study
%
% Options is the struct that contains the method options, the significant ones are:
%
% par.kvoi - Proportionality constant for the calculation of the
%            tracer concentration in the VOI, by default it is
%            considered unknown and set to 1.
%
% S0 - series of parameters to identify the S0 calculation threshold
% S0.nSamplesMin  - number of samples that are surely acquired
%                   before injection
% S0.nSamplesMax  - number of samples after which it stops in any case
% S0.thresh;      - I add a sample if its difference from the mean is
%                   less than the threshold
%
% options.display - level 1 Shows the progress of the processing
%                 - level 2 Shows the S0 maps and the average signal on
%                  which it was estimated
%
% Output parameters:
% - conc: 4D matrix of concentrations
% - S0map: 3D matrix of S0
% - figure 3: s0 maps and average signal on which it was estimated

function [conc,S0map,bolus]=DSC_mri_conc(volumes,mask,options)

if options.display > 0
    disp(' ')
    disp('Calculating concentration...');
end

[S0map,bolus]=DSC_mri_S0(volumes,mask,options);

conc=zeros(size(volumes));

ind=find(mask);
k=options.nR*options.nC*options.nS;
if options.waitbar
    hw=waitbar(0,'Calculating concentration...');
end




for t=1:options.nT
    if options.waitbar
        waitbar((t-1)/options.nT,hw)
    end
    step1=volumes(ind+k*(t-1))./S0map(ind); 
    %step 1 is a vector of x,1 where x is all pixels in the mask.
    % step1 is simply the raw signal in a pixel divided by the S0 value in that pixel.
    % the s0 value is the average signal intensity of the first few timepoints before the bolus has arrived.
    % so step1 is just the normalised signal intensity in a pixel.
 
    % disp(mean(step1(~isnan(step1))))
    % disp(t)

    % we would expect the step1 mean to be around 1, and then when the bolus arrives the step1 mean should be LOWER than 1.
    % This is because the raw signal has a dip when the bolus arrives.
    % In a good bolus the mean get to 0.8. So about a 20-25% drop in signal intensity.
    % In a bad bolus the mean is about 0.95. So about a 5% drop in signal intensity.
    % Then when log()=ln() in matlab. This function is only negative for values lower than 1.
    % we expect values lower than 1 (lower than the baseline). So there is also a negative sign in front of the scaling parameters.
    % negative and negative cancel out, so the conc values are positive. is the signal intensity is lower than the baseline, (Lower means more contrast agent.)

    % However when the bolus is weak the signal intensity is not lower than the baseline, combine this with slight patient movement and the conc value can jitter around 0 or get even lower. 
    % And low conc values lead, to a low integral of the conc curve. And a low integral of the conc curve when divided by a normal aif curve leads to a low cbv value.
    % And a cbv value around 0 leads to normalizing reference around zero which leads to exploding cbv values.
    
    %step2=step1.*(step1<1);
    %step3=step2+(step2==0);
    conc(ind+k*(t-1))=-(options.par.kvoi/options.te).*log(step1);
    % The concentration conc is calculated by normalizing the MRI signal at each time point with the baseline signal (S0map), and then applying a logarithmic transformation scaled by specific parameters (kvoi and te). The result is stored in the corresponding time slice of the conc matrix.
end

if options.waitbar
    delete(hw);
end



% Apply a Gaussian filter to the concentration maps
% This should smooth the concentration maps, which will hopefully reduce the noise in the signal per pixel due to jitter and slight patient motion.
% Use a 3x3 kernel in the slice plain. The axial resolution is much lower than the in-plane resolution so we dont want to smooth in the axial direction.

% disp('Applying Gaussian filter to concentration maps...');
% for s = 1:options.nS
%     for t = 1:options.nT
%         conc(:,:,s,t) = imgaussfilt(conc(:,:,s,t), 0.75);
%     end
% end

% % Calculate the number of negative values in the time (4th) dimension for each voxel
% negative_count = sum(conc < 0, 4);

% % Plot a middle slice of the negative count map
% figure;
% hAx = axes('Position', [0.1, 0.2, 0.8, 0.7]);
% hImg = imagesc(negative_count(:,:,round(0.5*options.nS)));
% colormap(gca, jet);
% colorbar;
% title('Number of negative concentration values per voxel');
% axis image;

% % Create slice scroll slider
% hSlider = uicontrol('Style', 'slider', ...
%     'Min', 1, ...
%     'Max', options.nS, ...
%     'Value', round(0.5*options.nS), ...
%     'SliderStep', [1/(options.nS-1), 1/(options.nS-1)], ...
%     'Units', 'normalized', ...
%     'Position', [0.1, 0.05, 0.8, 0.05]);

% % Add listener to update the image when the slider value changes
% addlistener(hSlider, 'Value', 'PostSet', @(src, event) updateNegativeCountImage(round(get(hSlider, 'Value')), negative_count, hImg));

% % Update function for the negative count image
% function updateNegativeCountImage(slice, negative_count, hImg)
%     set(hImg, 'CData', negative_count(:,:,slice));
%     title(hImg.Parent, ['Number of negative concentration values per voxel - Slice ' num2str(slice)]);
% end


end


%%
function [S0map,bolus]=DSC_mri_S0(volumes,mask,options)


% The function calculates the bolus injection time and S0 from the
% data.
% 1) on the average trend I calculate the bolus injection time: I calculate
%    the average of the first n samples and add the n+1 if its percentage
%    difference from the mean is less than a given threshold.
% 2) I calculate S0 as the average of the first n samples for all voxels.

% Define parameters
nSamplesMin=options.S0.nSamplesMin; % number of samples that are surely acquired before injection.
nSamplesMax=options.S0.nSamplesMax; % number of samples after which I stop in any case.
thresh=options.S0.thresh;

% --------------------------------------------------------------------------
% 1.1) calculate the average signal intensity for each slice and each time point
% --------------------------------------------------------------------------


% Create a matrix to store the mean signal for each slice and each time point.
% This matrix has size nSxnT. 
mean_signal=zeros(options.nS, options.nT);

if options.waitbar
    hw=waitbar(0,'Retrieving S0...');
end


for s=1:options.nS 
    for t=1:options.nT
        if options.waitbar
            waitbar((options.nT*(s-1)+(t-1))/(options.nT*options.nS),hw);
        end
        
        slice_data = volumes(:,:,s,t);
        %make the slice_data a vector
        slice_data = slice_data(:);
        % remove all the zeros
        slice_data = slice_data(slice_data~=0);
        % Get the mean of the slice_data
        mean_signal(s,t) = mean(slice_data);
    end
end


% A S0 map is calculated seperatly for each slice.
% This means that the injection time (start point) of the bolus is calculated separately for each slice. 
% The bolus injection is found by comparing the mean of that slice for an early  
% time point with the mean of the later time points.
% Once the difference between mean signal intensity between 2 consecutive time points is larger than the threshold,
% We assume that the bolus has arrived and we take that time point

% --------------------------------------------------------------------------
% 1.2) calculation of the bolus injection time
% --------------------------------------------------------------------------
% This is done individually for each slice.
for s=1:options.nS
    cycle=true;
    % pos is the time point which we are considering
    % pos is initialized to the minimum number of samples that are surely acquired before injection.
    % This saves a little time, since we know that the bolus has not arrived yet.
    pos=nSamplesMin;
    while cycle
        % Get the mean of the mean signal intensity up to the current time point "pos"
        mean_val=mean(mean_signal(s,1:pos));
        %if the difference between the mean up to the current time point and the mean up to the next time point is smaller than the threshold, 
        % we go on to the next time point.
        % If the difference exceeds the threshold, we assume that the bolus has arrived and we take that time point. 
        if abs((mean_val-mean_signal(s,pos+1))/mean_val)<thresh
            pos=pos+1;
        else
            cycle=false;
            pos=pos-1; % Conservative choice, do not consider
            % the last sample before injection
        end
        if pos==nSamplesMax
            cycle=false;
            pos=pos-1;
        end
    end

    % 2) Calculation of S0
    % each pixel of the s0 map is the mean signal intensity of the first few timepoints before the bolus has arrived.
    S0map(:,:,s)=mask(:,:,s).*mean(volumes(:,:,s,1:pos),4);
    bolus(s) = pos;

    
    
end

%--------------------------------------------------------------------------
% 1.3) Detect images with unclear bolus
%--------------------------------------------------------------------------
unclear_bolus = false(1, options.nS);
bolus_quality = cell(1, options.nS);

for s = 1:options.nS
    
    % as a first metric we calculate the signal drop from baseline to bolus dip
    % Get baseline (average signal before bolus) 
    baseline = mean(mean_signal(s, 1:bolus(s)));
    % Get the minimum (lowest point in the bolus dip)
    [min_signal, min_idx] = min(mean_signal(s, bolus(s):end));
    % Calculate the percentage signal drop from bolus dip to baseline
    signal_drop = (baseline - min_signal) / baseline * 100; 
    
    % As a second metric we calculate the signal to noise ratio.
    % The signal is absolute signal drop from baseline to bolus dip.
    % The noise is the standard deviation of the signal before bolus.
    % This snr shows if the bolus dip is significant compared to the random oscillations (noise)
    noise = std(mean_signal(s, 1:bolus(s)));
    snr = (baseline - min_signal) / noise;
    
    
    % Define thresholds
    min_drop_threshold = 10; %this is the minimum percentage drop in signal intensity that we consider a bolus, threshold was determined experimentally
    min_snr_threshold = 5;   % minimum snr that we count as a valid bolus.
 
    % if a slice fails one of the two quality metrics, we consider it unclear
    unclear_bolus(s) = (signal_drop < min_drop_threshold) || (snr < min_snr_threshold);
    
    % Store quality metric (higher is better)
    bolus_quality{s} = sprintf('%.2f%%, %.2f, %.2f', signal_drop, snr);

end

% Warning for unclear bolus, always print this.
if any(unclear_bolus) &&  options.display > 0
    disp(['WARNING: Unclear bolus detected in slices: ' num2str(find(unclear_bolus))]);
end

% if more than 50% of the slices are unclear, we print a warning to the error logs
if sum(unclear_bolus) > 0.5*options.nS
    % if "error_logs.txt" does not exist, create it
    if ~exist('error_logs.txt', 'file')
        fileID = fopen('error_logs.txt', 'w');
        fprintf(fileID, 'Error logs: \n');
        fclose(fileID);
    end
    
    % write the name of the file to the error_logs.txt file
    fileID = fopen('error_logs.txt', 'a');
    fprintf(fileID, 'Unclear bolus in at least 50 percent of slices for patient %s\n', options.perf_path);
    fclose(fileID);
end


if options.display >= 1
    disp(' ')
    disp('Bolus quality metrics per slice: (Percentage signal drop of bolus, SNR)');
    fprintf('Threshold for percentage signal drop is %.2f, threshold for SNR is %.2f\n', min_drop_threshold, min_snr_threshold);

    for s = 1:options.nS
        disp(['Slice ' num2str(s) ': ' bolus_quality{s}]);
    end
end


% Create plots after all slices are processed
if options.display > 1
    % Plot 1: Raw signal intensity map with time scroll
    middle_slice = round(0.5*options.nS);
    slice_w_timepoints = volumes(:,:,middle_slice,:);
    slice_w_timepoints = slice_w_timepoints .* mask(:,:,middle_slice);
    
    cmin = min(slice_w_timepoints(:));
    cmax = max(slice_w_timepoints(:));
    
    figure;
    hAx = axes('Position', [0.1, 0.2, 0.8, 0.7]);
    hImg = imagesc(slice_w_timepoints(:,:,1));
    colormap(gca, jet);
    caxis([cmin cmax]);
    title(['Raw signal for slice ' num2str(middle_slice)], [' Timepoint: 1']);
    colorbar;
    
    % Create time scroll slider
    hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(slice_w_timepoints, 4), ...
        'Value', 1, 'SliderStep', [1/(size(slice_w_timepoints, 4)-1), 1/(size(slice_w_timepoints, 4)-1)], ...
        'Units', 'normalized', 'Position', [0.1, 0.05, 0.8, 0.05]);
    
    addlistener(hSlider, 'Value', 'PostSet', @(src, event) updateImage(round(get(hSlider, 'Value')), slice_w_timepoints, hImg));
end

if options.display > 1	
    % Plot 2: S0 map and signal plot with slice scroll
    hf_s0 = figure;
    hAx1 = subplot(1, 3, 1); % Change to 1, 3 for three subplots
    
    % Calculate y-axis limits based on all slices data
    ymin = min(mean_signal(:));
    ymax = max(mean_signal(:));
    y_range = ymax - ymin;
    ylimits = [0, ymax + 0.05*y_range];
    
    hPlot = plot(options.time, mean_signal(middle_slice, :), 'b+-');
    hold on;
    hPoint = plot(options.time(bolus(middle_slice)), mean_signal(middle_slice, bolus(middle_slice)), 'ro');
    ylabel('Slice averaged raw signal intensity');
    xlabel('Time (s)');
    xlim([options.time(1) options.time(end)]);
    ylim(ylimits);  % Set consistent y-axis limits
    
    hAx2 = subplot(1, 3, 2); % Change to 1, 3 for three subplots
    hImg = imagesc(S0map(:, :, middle_slice));
    colormap(gca, jet);
    cmin = min(S0map(:));
    cmax = max(S0map(:));
    caxis([cmin cmax]);
    colorbar;
    axis image;
    xlabel('X-axis');
    ylabel('Y-axis');
    title({'S0 map', ['Slice ' num2str(middle_slice) '/' num2str(options.nS)], 'Slice averaged raw signal intensity'});
    impixelinfo;
    
    % New subplot for raw signal intensity at peak bolus point
    hAx3 = subplot(1, 3, 3); % New subplot
    peak_bolus_timepoint = bolus(middle_slice);
    hImgPeak = imagesc(volumes(:,:,middle_slice,peak_bolus_timepoint) .* mask(:,:,middle_slice));
    colormap(gca, jet);
    caxis([cmin cmax]);
    colorbar;
    axis image;
    xlabel('X-axis');
    ylabel('Y-axis');
    title({'Raw signal at peak bolus', ['Slice ' num2str(middle_slice) '/' num2str(options.nS)], ['Timepoint: ' num2str(peak_bolus_timepoint)]});
    impixelinfo;
    
    % Create slice scroll slider
    hSlider = uicontrol('Style', 'slider', ...
        'Min', 1, ...
        'Max', options.nS, ...
        'Value', middle_slice, ...
        'SliderStep', [1/(options.nS-1), 1/(options.nS-1)], ...
        'Units', 'normalized', ...
        'Position', [0.1, 0.05, 0.8, 0.05]);
    
    % Add listener with access to all necessary data
    addlistener(hSlider, 'Value', 'PostSet', @(src, event) updateSlice(min(max(round(get(hSlider, 'Value')), 1), options.nS), mean_signal, options, hPlot, hPoint, hImg, S0map, bolus, hImgPeak, volumes, mask));
end

if options.display > 1
    % plot the avearge signal intensity for the entire brain overtime, 

    mean_volumes = zeros(options.nT, 1);
    for time = 1:options.nT
        volume = volumes(:,:,:,time);
        volume_non_zero = volume(volume ~= 0);
        mean_volumes(time) = mean(volume_non_zero);
    end
 
    
    figure;
    plot(options.time, mean_volumes, 'b-');
    title('Average raw signal intensity in all brain voxels over time');
    xlabel('Time (s)');
    ylabel('Average signal intensity');
    legend('Average signal intensity');
    grid on;

end



if options.waitbar
    try
        delete(hw);
    end
end
end

% Add close function at end of file
function my_closefcn(src, ~)
    delete(src);
end

% Update helper functions
function updatePlotTitle(slice, pos, options, ax)
    title(ax, {['Slice ' num2str(slice) '/' num2str(options.nS)], 'Slice averaged raw signal intensity'});
end

function updateSlice(slice, mean_signal, options, hPlot, hPoint, hImg, S0map, bolus, hImgPeak, volumes, mask)
    % Ensure slice is within bounds
    slice = min(max(slice, 1), size(S0map, 3));
    
    % Update plot data
    set(hPlot, 'YData', mean_signal(slice, :));
    pos = bolus(slice);  % Use pre-calculated bolus position
    
    % Update point position
    set(hPoint, 'XData', options.time(pos), 'YData', mean_signal(slice, pos));
    
    % Keep y-axis limits consistent across slices
    ymin = min(mean_signal(:));
    ymax = max(mean_signal(:));
    y_range = ymax - ymin;
    ylim(hPlot.Parent, [0, ymax + 0.05*y_range]);
    
    % Update titles
    updatePlotTitle(slice, pos, options, hPlot.Parent);
    set(hImg, 'CData', S0map(:, :, slice));
    title(hImg.Parent, {'S0 map', ['Slice ' num2str(slice) '/' num2str(options.nS)], "S0 is the average until the bolus arrival (red dot)"});
    
    % Update peak bolus image
    [~, peak_bolus_timepoint] = min(mean_signal(slice, :));
    set(hImgPeak, 'CData', volumes(:,:,slice,peak_bolus_timepoint) .* mask(:,:,slice));
    title(hImgPeak.Parent, {'Raw signal at peak bolus', ['Slice ' num2str(slice) '/' num2str(options.nS)], ['Timepoint: ' num2str(peak_bolus_timepoint)]});
end

function updateImage(timepoint, slice_w_timepoints, hImg)
    set(hImg, 'CData', slice_w_timepoints(:,:,timepoint));
    title(hImg.Parent, ['Signal for slice, Timepoint: ' num2str(timepoint), "Slice averaged raw signal intensity"]);
end

