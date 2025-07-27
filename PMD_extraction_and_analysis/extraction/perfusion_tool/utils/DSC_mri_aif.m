function [aif]=DSC_mri_aif(conc,mask,options, atlas)
    % Input parameters:
    % - conc (4D Matrix) containing concentration time curves for all voxels
    % - mask (3D matrix) contains mask for each slice. Note: this mask is not
    %   the one used for concentration calculations, but a more restrictive version
    %   (fill function was not used)
    % - options is the struct containing method options, significant ones are:
    %
    % options.aif.enable: flag that enables AIF search (default=1)
    % 
    % options.aif.semiasseMaggiore: identifies AIF search area. The elliptical
    %                               search zone has major semi-axis equal to twice
    %                               the portion indicated by this option of the
    %                               brain size (default=0.35)
    %
    % options.aif.semiasseMinore: like major semi-axis, but for the other 
    %                            semi-axis of search region (default=0.15)
    %
    % options.aif.pArea: percentage of candidate voxels excluded based on
    %                   area under the curve (default=0.4)
    %
    % options.aif.pTTP: percentage of candidate voxels excluded based on
    %                  time to peak (default=0.4)
    %
    % options.aif.pReg: percentage of candidate voxels excluded based on
    %                  curve irregularity (default=0.05)
    %
    % options.aif.diffPicco: 0.0400
    %
    % options.aif.nVoxelMax: maximum number of accepted arterial voxels
    %                       (default=6)
    %
    % options.aif.nVoxelMin: minimum number of accepted arterial voxels  
    %                       (default=4)
    %
    % Output parameters: 
    % - aif structure, containing ROI
    % - figure 4: left: the sliced slice for AIF extraction, with the selected voxels used to calculate the AIF
    % - figure 4: right: the signal intensity for the selected voxels and the GV (Gamma variate) fit that averages the signal intensities 
    
    
    %%%%%%%%%%%%%%%%%%% MAIN CODE TO CALL A FUNCTION THAT EXTRACT THE AIF FROM EVERY SLICE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    addpath('utils_aif');
    
    if options.display > 0
        disp(' ')
        disp('AIF extraction...');
    end

    % Save a deep copy of the original concentration data
    conc_original = conc;
    
    % Initialse a list to store the aif results for each slice
    aif_plot_list = cell(options.nS,1);
    
    % Make a mask of the reference region in the atlas
    % With this atlas you make sure the aif can only be extracted for a specific brain region to make it more concistent across patients.
    atlas_ref_region = zeros(size(atlas));
    atlas_ref_region(atlas == 29 | atlas == 30 | atlas == 77 | atlas == 78 ) = 1;
    % These 4 regions make up the cingulate gyrus. This regions was chosen because it is fairly central in the brain and large vessels of the ACA run through it.

    % for every slice, extract the aif results and add them to aif_plot_list
    for slice_id = 1:options.nS
        
        %print a whiteline to separate the slices better when printing the output

        if options.display > 0	
            disp(' ')
        end


        reshaped_conc_slice = conc_original(:,:,slice_id,:);
        reshaped_conc_slice = squeeze(reshaped_conc_slice); %remove the singleton dimension
        mask_slice = mask(:,:,slice_id);
        atlas_ref_slice = atlas_ref_region(:,:,slice_id);
    
        % Check if the mask slice is empty, in that cause you cannot extract the AIF, so we add dummy data to the lists, such that the correspondence to the number of slices is maintained.
        % the same goes for slices where the mask is very small, in the top few slices, we will not find the suitable aif there either,
        % So we set threshold for the amount of voxels in the mask, this should be based on the specific dataset. but i would make it similar to the roi size in number of voxels.
        if sum(mask_slice(:)) < options.min_mask_size
            % Skip extraction and add dummy data to the lists
            aif_slice_plot_data = struct('image', [], 'bound', [], 'ROI_x', [], 'ROI_y', [], 'voxels', [], 'time', [], 'conc', [], 'fit_gv', []);
            aif_old = struct('conc', [], 'fit', struct('gv', []));
            if options.display > 0
                disp('slice_id: ' + string(slice_id) + ' brain mask too small, skipping AIF extraction');
            end
        else
            % Proceed with AIF extraction
            [aif_old, aif_slice_plot_data] = extractAIF(reshaped_conc_slice, mask_slice, options, atlas_ref_slice);
            if options.display > 0
                disp(['AIF extraction completed for slice ' num2str(slice_id)]);	
            end
        end
                
        aif_plot_list{slice_id} = aif_slice_plot_data;
        aif_old_list{slice_id} = aif_old;
    end
    
   

    %%%%%%%%%%%%%%%%%%%%%%%% CODE TO GENERATE A FIGURE WHERE YOU CAN SCROLL THROUGH THE AIF RESULTS FOR EVERY SLICE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.display > 1
        % Create a figure with a slider to navigate through slices
        hf_slider = figure('Name', 'AIF Extraction with Slider');
    
        % Create axes for anatomical view and AIF curves
        ax1 = subplot(1, 2, 1);
        ax2 = subplot(1, 2, 2);
    
        % Create slider
        slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', options.nS, 'Value', 1, ...
                        'Units', 'normalized', 'Position', [0.3 0.01 0.4 0.05], ...
                        'SliderStep', [1/(options.nS-1) 1/(options.nS-1)], ...
                        'Callback', @updateSlice);
    
        % Create text to display current slice number
        slice_text = uicontrol('Style', 'text', 'Units', 'normalized', ...
                            'Position', [0.75 0.01 0.1 0.05], ...
                            'String', ['Slice: 1/' num2str(options.nS)]);
    
        % Initial plot
        updateSlice(slider);
    end
    
    % Callback function to update the slice
    function updateSlice(hObject, ~)
        slice_id = round(get(hObject, 'Value'));
        set(slice_text, 'String', ['Slice: ' num2str(slice_id) '/' num2str(options.nS)]);
        
        if isempty(aif_plot_list{slice_id}.image)
            % Skip plotting if the image is empty
            cla(ax1);
            cla(ax2);
            return;
        end
        
        % Update anatomical view
        axes(ax1);
        imagesc(real(aif_plot_list{slice_id}.image), aif_plot_list{slice_id}.bound), colormap(gray)
        hold on
        plot(aif_plot_list{slice_id}.ROI_x, aif_plot_list{slice_id}.ROI_y, 'r', 'LineWidth', 1)
        plot(aif_plot_list{slice_id}.voxels(:,2), aif_plot_list{slice_id}.voxels(:,1), 'r.', 'MarkerSize', 2)
        xlabel('', 'FontSize', 10)
        legend({'Searching area', 'Selected voxels'}) % Ensure legend entries are in a cell array
        title(['Arterial voxel location - Slice ' num2str(slice_id)], 'FontSize', 12)
        set(gca, 'xtick', [], 'ytick', [], 'fontsize', 10)
        axis square
        hold off
        
        % Update AIF curves
        axes(ax2);
        plot(aif_plot_list{slice_id}.time, aif_plot_list{slice_id}.conc, 'ko', 'MarkerSize', 5)
        hold on
        plot(aif_plot_list{slice_id}.time, aif_plot_list{slice_id}.fit_gv, 'k-', 'LineWidth', 2)
        xlabel('time [s]', 'FontSize', 10)
        ylabel('Concentration ([mmol/L]?)', 'FontSize', 10)
        
        ylim([-5 100]) %preferably the 100 here is changed to be the max data for the aif across all slices time 1.25 for some extra room.
        xlim([0 size(conc, 4)]) 
        legend({'AIF samples', 'Gamma variate fit'}) % Ensure legend entries are in a cell array
        title(['AIF - Slice ' num2str(slice_id)], 'FontSize', 12)
        xlim([aif_plot_list{slice_id}.time(1) aif_plot_list{slice_id}.time(end)])
        hold off
    end



%%%%%%%%%%%%%%%%%% SELECT THE SLICE WITH THE BEST AIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create a data structure to store peak data and MSE for each slice
    slice_ids = num2cell(1:options.nS)';            
    peak_heights = num2cell(zeros(options.nS, 1));  
    mses = num2cell(zeros(options.nS, 1));          
    
    
    peak_data = struct('slice_id', slice_ids, ...
                      'peak_height', peak_heights, ...
                      'mse', mses);
    
    % Populate the data structure with peak data for each slice
    for slice_id = 1:options.nS
        if ~isempty(aif_plot_list{slice_id}.fit_gv)
            
            % Find the peak height 
            peak_data(slice_id).peak_height = max(aif_plot_list{slice_id}.fit_gv);
            
            % Calculate the mean square error between the fit and the data
            peak_data(slice_id).mse = mean((aif_plot_list{slice_id}.conc - aif_plot_list{slice_id}.fit_gv).^2);
        else
            peak_data(slice_id).peak_height = -Inf; % Assign a very low value if fit_gv is empty
            peak_data(slice_id).mse = Inf; % Assign a very high value if fit_gv is empty
        end
    end
    

    % Sort the data structure based on peak height
    [~, sorted_peak_height_indices] = sort([peak_data.peak_height], 'descend');
    sorted_peak_height_data = peak_data(sorted_peak_height_indices);
    
    % Sort the data structure based on mse
    [~, sorted_mse_indices] = sort([peak_data.mse]);
    sorted_mse_data = peak_data(sorted_mse_indices);
    
    
    if options.display > 0
        % Display the sorted lists
        disp('Slices sorted by peak height:');
        disp(' ')
        for i = 1:options.nS
            disp(['Slice ' num2str(sorted_peak_height_data(i).slice_id) ': Peak Height = ' num2str(sorted_peak_height_data(i).peak_height) ', mse = ' num2str(sorted_peak_height_data(i).mse)]);
        end
        disp('From the top three slices based on peak height, the best slice is selected based on the lowest mse.');
    end
    

    % Take only the 3 best slices based on the peak height
    best_slices_based_height = [sorted_peak_height_data(1:3)];

    % Sort the best slices based on mse
    [~, best_mse_indices] = sort([best_slices_based_height.mse]);
    sorted_peak_then_mse = best_slices_based_height(best_mse_indices);

    % Take the best slice based on the mse
    best_slice_id = sorted_peak_then_mse(1).slice_id;
    if options.display > 0
        disp(['Best slice selected: ' num2str(best_slice_id)]);
    end
    
    % Extract AIF from the best slice
    aif = aif_old_list{best_slice_id};
    aif.sliceAIF = best_slice_id;

    % Apply quadratic relaxivity correction if enabled
    if options.qr.enable
        % Convert AIF concentration using quadratic relationship:
        % C_new = (a*C + b*C^2)/r 
        % where a,b,r are relaxivity parameters
        aif.conc = (options.qr.a*aif_old.conc + ...
                    options.qr.b*(aif_old.conc.^2))/options.qr.r;
        
        % Display comparison plots if verbose mode enabled            
        if options.display > 1
            hf.aif = figure();
            disp('Converting AIF concentration using quadratic relationship ')
    
            subplot(121)
            plot(options.time, aif_old.conc)
            title('AIF with linear \Delta R_2^* relationship')
            
            subplot(122)
            plot(options.time, aif.conc) 
            title('AIF with quadratic \Delta R_2^* relationship')
            
            % Wait for the user to close the figure
            set(hf.aif, 'CloseRequestFcn', @my_closefcn);
        end
    end
    
    
    function [AIF, aif_slice_plot_data]=extractAIF(AIFslice,mask,options, atlas_ref_slice)
    % Extracts the AIF from the provided slice
    %
    % This function identifies the  Arteral Input Function (AIF) from a given MRI slice. 
    % The methodology involves several steps to ensure accurate selection of arterial voxels. 
    % The process begins by defining a Region of Interest (ROI) using an elliptical mask based on brain dimensions. 
    % The major and minor semi-axes of the ellipse are determined by user-defined parameters. 
    % The ROI is further refined by excluding voxels outside the brain mask.
    %
    % The next step involves decimating (removing) candidate voxels based on three criteria: 
    % Area Under the Curve (AUC), Time To Peak (TTP), and curve irregularity. 
    % Initially, voxels are filtered based on their AUC, retaining those with higher concentration values. 
    % A binary search algorithm is used to find an AUC threshold that keeps a specified percentage of voxels. 
    % The remaining voxels are then filtered based on their TTP, with a threshold determined similarly to the AUC step. 
    % Finally, voxels are filtered based on their irregularity index, calculated from the second derivative of the concentration curves. 
    % This step ensures that smoother curves, indicative of true arterial signals, are retained.
    %
    % After decimation, a hierarchical clustering algorithm is applied to the remaining voxels. 
    % The algorithm divides the voxels into two clusters, and the cluster with the higher peak concentration or earlier TTP 
    % is selected as the arterial cluster. This process is repeated recursively until the number of voxels falls below a user-defined maximum.
    %
    % The final step involves preparing the output, which includes the indices of the selected voxels, their coordinates, and the mean concentration curve. 
    % A gamma-variate fit is applied to the concentration curve to model the AIF, with optional fitting for recirculation. 
    % The fitting parameters and the fitted curve are included in the output structure.
    %
    % The function also includes options for visualizing intermediate and final results, aiding in the verification of the AIF extraction process.
    
    
    % Preparation of auxiliary variables and parameters
    semiasseMag  = options.aif.semiasseMaggiore; % Major semi-axis for search ellipse 
    semiasseMin  = options.aif.semiasseMinore;   % Minor semi-axis for search ellipse
    pArea = options.aif.pArea;         % % of voxels to exclude based on area
    pTTP = options.aif.pTTP;          % % of voxels to exclude based on time-to-peak
    pReg = options.aif.pReg;          % % of voxels to exclude based on irregularity
    nVoxelMax = options.aif.nVoxelMax; % Max number of arterial voxels to accept
    nVoxelMin = options.aif.nVoxelMin; % Min number of arterial voxels to accept
    diffPicco = options.aif.diffPicco;  % Allowed difference between peak heights
    
    % Get dimensions of input slice
    [nR,nC,nT]=size(AIFslice);
    maskData=mask;    % Store mask in separate variable
    
    % Prepare image for visualizations
    % Sum up slice over time dimension for display
    image.img=sum(AIFslice,3);
    % Calculate display bounds from 95th percentile
    vectorImage=sort(image.img(1:nR*nC)); 
    image.bound=[0 vectorImage(round(0.95*nR*nC))];
    %sometimes is the area in the mask is very small, like for the last slices, the 95 percentile range gives bounds 0 0
    % In that case use the full range of the image
    if image.bound(2) == 0
        image.bound(2) = max(vectorImage);
    end
    clear vectorImage
    
    % --------------------------------------------------------------
    % 1) IDENTIFICATION OF ROI CONTAINING THE AIF
    % --------------------------------------------------------------
    % 1.1) Find mask boundaries
    
    if options.display > 2
        disp('   Brain bound detection')
    end
    
    % Find minimum row with signal
    cycleFlag=true;
    r=1;
    while cycleFlag
        if sum(maskData(r,:))~=0
            minR=r;
            cycleFlag=false;
        else
            r=r+1;
        end
    end
    
    % Find maximum row with signal
    cycleFlag=true; 
    r=nR;
    while cycleFlag
        if sum(maskData(r,:))~=0
            maxR=r;
            cycleFlag=false;
        else
            r=r-1;
        end
    end
    
    % Find minimum column with signal
    cycleFlag=true;
    c=1;
    while cycleFlag
        if sum(maskData(:,c))~=0
            minC=c;
            cycleFlag=false;
        else
            c=c+1;
        end
    end
    
    % Find maximum column with signal
    cycleFlag=true;
    c=nC;
    while cycleFlag
        if sum(maskData(:,c))~=0
            maxC=c;
            cycleFlag=false;
        else
            c=c-1;
        end
    end
    
    % Display mask boundaries if verbose mode enabled
    if options.display > 2
        hf.mask=figure();
        imagesc(maskData)
        hold on
        % Plot boundary lines
        plot([1 options.nC],[minR minR],'g-',[1 options.nC],[maxR maxR],'g-',...
             [minC minC],[1 options.nR],'g-',[maxC maxC],[1 options.nR],'g-')
        xlabel(['Brain bound: rows (' num2str(minR) '-' num2str(maxR) ...
                ') - columns(' num2str(minC) '-' num2str(maxC) ')'])
        title('AIF extraction - mask and bounds')
        set(gca,'xtick',[],'ytick',[])
        axis square
        set(hf.mask, 'CloseRequestFcn', @my_closefcn);
    end
    
    
    
    
    
    % 1.2) Drawing of the Region of Interest (ROI)
    if options.display > 2
        disp('   Definition of the AIF extraction searching area')
    end
    
    % Calculate center coordinates of brain mask
    centro(2)=0.5*(minR+maxR); % Y coordinate of center (calculated from rows)
    centro(1)=0.5*(minC+maxC); % X coordinate of center (calculated from columns)
    
    % The major semi-axis is along the anterior-posterior direction (left to right in images)
    % Calculate size of elliptical search region based on brain dimensions
    semiasseB=semiasseMag.*(maxC-minC); % Major semi-axis length
    semiasseA=semiasseMin.*(maxR-minR); % Minor semi-axis length
    
    % Initialize ROI mask
    ROI=zeros(nR,nC); % Mask containing ROI voxels
    
    % Create elliptical ROI using equation (x-h)^2/a^2 + (y-k)^2/b^2 <= 1
    % where (h,k) is center, a is semi-major axis, b is semi-minor axis 
    for r=1:nR
        for c=1:nC
            if ((r-centro(2))^2/(semiasseA^2)+(c-centro(1))^2/(semiasseB^2))<=1
                ROI(r,c)=1;
            end
        end
    end
    
    % Keep only voxels that are in both ROI and brain mask, AND IN THE ATLAS REF REGION
    ROI=ROI.*maskData.*atlas_ref_slice; 
    % ROI=ROI.*maskData; 
    ROIiniziale=ROI; % Store initial ROI for later use
    
    % Generate points to draw ellipse outline for visualization
    % Create x coordinates with fine sampling
    xROI=centro(1)-semiasseB:0.01:centro(1)+semiasseB;
    nL=length(xROI);
    xROI(2*nL)=0;
    yROI=zeros(1,2*nL);
    
    % Calculate upper half of ellipse
    for k=1:nL
        yROI(k)=semiasseA*((1-((xROI(k)-centro(1))^2)/(semiasseB^2))^0.5)+centro(2);
    end
    
    % Calculate lower half of ellipse by mirroring upper half
    for k=1:nL
        xROI(nL+k)=xROI(nL-k+1);
        yROI(nL+k)=-semiasseA*((1-((xROI(nL+k)-centro(1))^2)/(semiasseB^2))^0.5)+centro(2);
    end
    
    % Remove any imaginary components from numerical errors
    yROI=real(yROI);
    
    % Display ROI overlay on brain image if verbose mode enabled
    if options.display > 2
        hf.img_roi=figure();
        imagesc(image.img,image.bound),colormap(gray)
        hold on
        plot(xROI,yROI,'r') % Draw ellipse outline in red
        plot(centro(1),centro(2),'r+') % Mark center point
        title('AIF extraction - searching area')
    end
    
    
    %--------------------------------------------------------------
    % 2) DECIMATION OF CANDIDATE VOXELS
    %--------------------------------------------------------------
    if options.display > 2
        disp('   Candidate voxel analysis')
    end
    
    % 2.1) AUC based selection of voxels which could potentially be used to calcuate the AIF.
    % Selection of suitable voxels is based on area under the signal curve (AUC) of the 4D concentration matrix.
    % This part filters voxels based on their total AUC to keep only those
    % with higher concentration values that are more likely to be arterial.
    
    % Calculate how many voxels to keep based on pArea parameter 
    totCandidates = sum(sum(ROI));
    totCandidatesToKeep = ceil(totCandidates.*(1-pArea));
    
    % Calculate AUC for each voxel by summing over time dimension
    AUC = sum(AIFslice,3); 
    AUC = AUC.*ROI; % Only keep AUC values within ROI
    AUC(isinf(AUC)) = 0; % Remove any infinite values
    
    % Binary search to find AUC threshold that keeps desired number of voxels
    cycleFlag = true;
    nCycle = 0;
    AUCdown = min(min(AUC)); % Lower bound
    AUCup = max(max(AUC)); % Upper bound
    
    while cycleFlag
        nCycle = nCycle + 1;
        threshold_auc = 0.5*(AUCup+AUCdown); % Test threshold halfway between bounds
        nCandidates = sum(sum(AUC>threshold_auc)); % Count voxels above threshold
        
        % Update bounds based on whether we have too many or too few voxels
        if nCandidates == totCandidatesToKeep
            cycleFlag = false;
        elseif nCandidates > totCandidatesToKeep
            AUCdown = threshold_auc; % Need higher threshold
        else
            AUCup = threshold_auc; % Need lower threshold
        end
        
        % Exit if bounds are very close or too many iterations
        if ((AUCup-AUCdown)<0.01)||(nCycle>100)
            cycleFlag = false;
        end
    end
    
    % Create mask showing kept vs rejected voxels:
    % 2 = rejected voxel
    % 1 = kept voxel
    ROIauc = 2.*ROI - ROI.*(AUC>threshold_auc); 
    
    % Update ROI to only keep voxels above AUC threshold
    ROI = ROI.*(AUC>threshold_auc);
    
    
    % Display results if requested
    if options.display > 2 
        disp(' ')
        disp(' Candidate voxel selection via AUC criteria')
        disp(['  Initial voxel count: ' num2str(totCandidates)])
        disp(['  Surviving voxels: ' num2str(sum(sum(ROI)))])
        
        % Create figure showing:
        % Top left: Spatial map of kept/rejected voxels
        % Bottom left: Time curves for all voxels (red=kept, blue=rejected)
        hf.aif = figure();
        subplot(2,3,1)
        imagesc(ROIauc)
        title('AUC criteria - survived voxels')
        set(gca,'xtick',[],'ytick',[])
        axis square
        
        subplot(2,3,4)
        plot(options.time,zeros(1,nT),'b-')
        hold on
        plot(options.time,zeros(1,nT),'r-')
        plot(options.time,zeros(1,nT),'k-')
        
        % Plot rejected voxel curves in blue
        for c=1:nC
            for r=1:nR
                if ROIauc(r,c)==2
                    plot(options.time,reshape(AIFslice(r,c,:),1,nT),'b-')
                end
            end
        end
        
        % Plot kept voxel curves in red
        for c=1:nC
            for r=1:nR
                if ROIauc(r,c)==1
                    plot(options.time,reshape(AIFslice(r,c,:),1,nT),'r-')
                end
            end
        end
        
        set(gca,'FontSize',12)
        legend('Accepted','Rejected')
        title('AUC')
        xlabel('time')
    end
    
    
    % 2.2) Selection based on Time To Peak (TTP)
    % This section filters voxels based on how quickly they reach their peak concentration
    
    % Calculate how many voxels to keep
    totalCandidates = sum(sum(ROI));
    totalCandidatesToKeep = ceil(totalCandidates.*(1-pTTP));
    
    % Find time to peak for each voxel by finding max along time dimension
    [MC,TTP] = max(AIFslice,[],3); 
    TTP = TTP.*ROI; % Only keep TTP values within ROI
    
    % Binary search to find TTP threshold that keeps desired number of voxels
    keepSearching = true;
    threshold_ttp = 1;
    while keepSearching
        % Count how many voxels have TTP less than threshold (excluding zeros)
        if (sum(sum(TTP<threshold_ttp))-sum(sum(TTP==0)))>=totalCandidatesToKeep
            keepSearching = false;
        else
            threshold_ttp = threshold_ttp+1;
        end
    end
    
    % Create mask: 2 = rejected voxel, 1 = kept voxel
    ROIttp = 2*ROI-ROI.*(TTP<threshold_ttp); 
    
    % Update ROI to only keep voxels below TTP threshold
    ROI = ROI.*(TTP<threshold_ttp);
    
    
    % Display results if verbose mode enabled
    if options.display > 2
        disp(' ')
        disp(' Candidate voxel selection via TTP criteria')
        disp(['  Initial voxel count: ' num2str(totalCandidates)])
        disp(['  Surviving voxels: ' num2str(sum(sum(ROI)))])
        
        % Plot spatial map of kept/rejected voxels
        figure(hf.aif);
        subplot(232)
        imagesc(ROIttp)
        title('TTP criteria - survived voxels')
        set(gca,'xtick',[],'ytick',[])
        axis square
        
        % Plot time curves for kept and rejected voxels
        subplot(235)
        plot(options.time,zeros(1,nT),'b-')
        hold on
        plot(options.time,zeros(1,nT),'r-')
        plot(options.time,zeros(1,nT),'k-')
        for c=1:nC
            for r=1:nR
                % Plot kept voxels in red, rejected in blue
                if ROIttp(r,c)==1
                    plot(options.time,reshape(AIFslice(r,c,:),1,nT),'r-')
                elseif ROIttp(r,c)==2
                    plot(options.time,reshape(AIFslice(r,c,:),1,nT),'b-')
                end    
            end    
        end    
        set(gca,'FontSize',12)
        legend('Accepted','Rejected')
        title('TTP')
        xlabel('time')
    end    
    
    
    
    
    % 2.3) Selection based on curve irregularity index
    % This section filters voxels based on how irregular their concentration curves are
    % Smoother curves are preferred as they likely represent true arterial signals
    
    totalCandidates = sum(sum(ROI));
    totalCandidatesToKeep = ceil(totalCandidates.*(1-pReg));
    
    % Use the calcolaReg function from the utils file
    REG = calcolaReg(AIFslice,options.time,ROI);
    
    % Binary search to find irregularity threshold
    keepSearching = true;
    nCycle = 0;
    REGdown = min(min(REG));
    REGup = max(max(REG));
    
    while keepSearching
        nCycle = nCycle+1;
        threshold_reg = 0.5*(REGup+REGdown);
        nCandidates = sum(sum(REG>threshold_reg));
        
        if nCandidates==totalCandidatesToKeep
            keepSearching = false;
        elseif nCandidates<totalCandidatesToKeep
            REGup = threshold_reg;
        else
            REGdown = threshold_reg;
        end
        
        % Exit if bounds are very close or too many iterations
        if ((REGup-REGdown)<0.001)||(nCycle>=100)
            keepSearching = false;
        end
    end
    
    % Create mask: 2 = rejected voxel, 1 = kept voxel
    ROIreg = 2*ROI-ROI.*(REG>threshold_reg);
    
    % Update ROI to only keep voxels above irregularity threshold
    ROI = ROI.*(REG>threshold_reg);
    

    % Call the function to display results if verbose mode enabled
    if options.display > 2
        displayVerboseResults(options, totalCandidates, ROI, ROIreg, hf, AIFslice, image, xROI, centro, ROIiniziale, AUC, threshold_auc, TTP, threshold_ttp, REG, threshold_reg, nT, nC, nR);
    end
    
    % --------------------------------------------------------------
    % 3) APPLICATION OF CLUSTERING ALGORITHM TO FIND ARTERIAL VOXELS
    % --------------------------------------------------------------
    if options.display > 2
        disp('   Arterial voxels extraction')
    end
    
    % 3.1) Prepare matrix containing the data
    % Convert 4D data to 2D matrix where each row represents a voxel's time curve
    dati2D = zeros(sum(sum(ROI)),nT); % Initialize matrix
    ind = find(ROI); % Find indices of voxels in ROI
    
    
    % Check if there are any voxels in the ROI
    if options.display > 0
        disp(['Number of voxels in the ROI of the atlas: ' num2str(length(ind))]);
    end
    if length(ind) < 2
        aif_slice_plot_data = struct('image', [], 'bound', [], 'ROI_x', [], 'ROI_y', [], 'voxels', [], 'time', [], 'conc', [], 'fit_gv', []);
        aif_old = struct('conc', [], 'fit', struct('gv', []));
        AIF = 0;
        if options.display > 0
            disp('We need at least 2 aterial voxels in the atlas ROI to perform clustering, now skipping this slice');
        end
        return;
    end
    
    
    
    k = nR*nC;
    for t=1:nT
        dati2D(:,t) = AIFslice(ind+k*(t-1)); 
    end
    maskAIF = ROI; % Initialize mask for arterial voxels
    
    % 3.2) Apply Hierarchical Clustering algorithm recursively
    % Keep dividing clusters until we get small enough group of arterial voxels
    
    cycleFlag = true;
    nCycle = 0;
    clear AIFslice
    
    while cycleFlag
        nCycle = nCycle + 1;
        if options.display > 2
            disp(' ')
            disp(' ------------------------------------')
            disp(['  CYCLE N# ' num2str(nCycle)])
        end
        
        % Use the clusterGerarchico function from the utils file
        [clusterAssignments, centroids] = clusterGerarchico(dati2D,2);
        
        % Compare the clusters and choose which one to keep
        % Find peak height and time-to-peak for each cluster
        [peakHeight1,TTP1] = max(centroids(1,:)); % Cluster 1 
        [peakHeight2,TTP2] = max(centroids(2,:)); % Cluster 2
        
        % Decide which cluster represents arterial voxels based on:
        % 1) If peak heights are similar (within diffPicco threshold), choose based on earlier TTP
        % 2) Otherwise choose cluster with higher peak
        if (((max([peakHeight1 peakHeight2])-min([peakHeight1 peakHeight2]))/max([peakHeight1 peakHeight2]))<diffPicco)&&(TTP1~=TTP2)
            % Peak difference is smaller than threshold, select based on time-to-peak
            chosenCluster = 1+(TTP2<TTP1); % Choose cluster 1 if TTP1<TTP2, cluster 2 otherwise
            
            if options.display > 2
                disp('  Cluster selected via TTP criteria')
                disp(['   Selected cluster: ' num2str(chosenCluster)])
            end
            
        else
            % Select based on peak height difference
            chosenCluster = 1+(peakHeight2>peakHeight1); % Choose cluster 1 if MC1>MC2, cluster 2 otherwise
            
            if options.display > 2
                disp('  Cluster selected via peak height criteria')
                disp(['   Selected cluster: ' num2str(chosenCluster)])
            end
        end
        
        % Check if chosen cluster has enough voxels
        % If chosen cluster has fewer than minimum required voxels but other cluster has enough,
        % switch to the other cluster
        if (sum(clusterAssignments==chosenCluster)<nVoxelMin)&&(sum(clusterAssignments==(3-chosenCluster))>=nVoxelMin)
            chosenCluster = 3-chosenCluster; % Switch clusters (1->2 or 2->1)
            
            if options.display > 2
                disp('  Cluster selection switched due to minimum voxel requirement')
                disp(['   Selected cluster: ' num2str(chosenCluster)])
            end
        end
        
        % Keep only the data from chosen cluster
        selectedVoxels = (clusterAssignments==chosenCluster);
        maskIndices = find(maskAIF);
        maskAIF(maskIndices) = selectedVoxels;
        
        % Update data matrix to only include chosen voxels
        voxelIndices = find(selectedVoxels);
        nVoxels = length(voxelIndices);
        dati2Dold = dati2D;
        dati2D = zeros(nVoxels,nT);
        for t=1:nT
            dati2D(:,t) = dati2Dold(voxelIndices,t);
        end
        
        % Display progress if verbose mode enabled
        if options.display > 2
            disp(' ')
            disp([' Resume cycle n# ' num2str(nCycle)])
            disp(['  Initial voxel count: ' num2str(length(maskIndices))])
            disp(['  Surviving voxels: ' num2str(nVoxels)])
            disp(['  Cluster 1: Peak height ' num2str(peakHeight1)])
            disp(['             TTP      ' num2str(TTP1)])
            disp(['             voxels    ' num2str(sum(clusterAssignments==1))])
            disp(['  Cluster 2: Peak height ' num2str(peakHeight2)])
            disp(['             TTP      ' num2str(TTP2)])
            disp(['             voxels    ' num2str(sum(clusterAssignments==2))])
            disp(['  Selected cluster: ' num2str(chosenCluster)])
            
            % Create visualization of current clustering results
            eval(['hf.img_centr' num2str(nCycle) '=figure();']);
            subplot(1,2,1)
            [posC,posR] = find(maskAIF);    
            imagesc(image.img,image.bound)
            colormap(gray)
            hold on
            plot(xROI,yROI,'r')
            plot(posR,posC,'r.','MarkerSize',1)
            title(['Cycle n#' num2str(nCycle) ' - candidate voxels'])
            set(gca,'xtick',[],'ytick',[],'fontsize',12)
            axis square
            
            subplot(1,2,2)
            plot(options.time,centroids,'k-')
            hold on
            plot(options.time,centroids(chosenCluster,:),'r-')
            title('Cluster centroids')
            xlabel('time')
            set(gca,'fontsize',12)
            eval(['set(hf.img_centr' num2str(nCycle) ', ''CloseRequestFcn'', @my_closefcn);']);
        end
        
        % Check exit criteria:
        % Stop if we have few enough voxels or hit maximum iterations
        if (nVoxels<=nVoxelMax)||(nCycle>=100)
            cycleFlag = false;
        end
    end
    
    % --------------------------------------------------------------
    % 4) PREPARATION OF OUTPUT
    % --------------------------------------------------------------
    
    % 4.1) Save the search ROI information
    % Store indices of voxels within ROI and the ellipse boundary coordinates
    AIF.ROI.ind = find(ROI); 
    AIF.ROI.x = xROI; % x-coordinates of elliptical search region
    AIF.ROI.y = yROI; % y-coordinates of elliptical search region
    
    % 4.2) Save position of chosen voxels and mean concentration
    % Get concentration time curves from selected cluster centroid 
    AIFconc = centroids(chosenCluster,:); % AIF concentration samples
    AIF.conc = AIFconc;
    
    % Store row/column coordinates of selected arterial voxels
    pos = 1;
    for r=1:nR
        for c=1:nC
            if maskAIF(r,c)==1
                AIF.voxels(pos,1) = r;
                AIF.voxels(pos,2) = c;
                pos = pos+1;
            end
        end
    end
    
    % 4.3) Calculate gamma-variate fit of arterial curve (with recirculation)
    if options.display > 2
        disp('   Gamma variate fit computation')
    end
    
    % Calculate weights for fitting
    % Lower weights for high concentration values to reduce their influence
    pesi = 0.01 + exp(-AIFconc); 
    
    % Adjust weights around the peak to improve fit quality
    [MC, TTP] = max(AIFconc);
    % if the ttp is below 3, it will gives us some errors because later the code assumes 
    % that ttp is atleast 3. So if the ttp for a slice is less than 3, 
    % we will skip that slice and return empty values.
    if TTP < 3
        aif_slice_plot_data = struct('image', [], 'bound', [], 'ROI_x', [], 'ROI_y', [], 'voxels', [], 'time', [], 'conc', [], 'fit_gv', []);
        aif_old = struct('conc', [], 'fit', struct('gv', []));
        disp('TTP is less than 3, this is probably an error, so we will skip this slice');
        return;
    end
    
    
    % Check if TTP+1 is within the range of pesi
    if TTP + 1 > length(pesi)
        aif_slice_plot_data = struct('image', [], 'bound', [], 'ROI_x', [], 'ROI_y', [], 'voxels', [], 'time', [], 'conc', [], 'fit_gv', []);
        aif_old = struct('conc', [], 'fit', struct('gv', []));
        disp('Warning: TTP+1 exceeds the range of pesi.this is probably an error, so we will skip this slice');
        return;
    end
    
    pesi(TTP) = pesi(TTP)./10;     % Reduce weight at peak
    pesi(TTP-1) = pesi(TTP-1)./5;  % Reduce weight before peak
    pesi(TTP+1) = pesi(TTP+1)./2;  % Reduce weight after peak
    
    % Parameter descriptions for output
    p = {'t0' ''; 'alpha' ''; 'beta' ''; 'A' ''; 'td' ''; 'K' ''; 'tao' ''; 'ExitFlag' ''};
    
    % Fit first pass (main bolus peak)
    [fitParameters_picco1,cv_est_parGV_picco1] = fitGV_picco1(AIFconc,pesi,options);
    
    % If recirculation fitting is enabled, fit second pass
    if options.aif.ricircolo
        % Fit recirculation peak using residuals from first pass
        [fitParameters_picco2,cv_est_parGV_picco2] = fitGV_picco2(AIFconc,pesi,fitParameters_picco1,options);
        
        % Combine parameters from both fits
        fitParameters = [fitParameters_picco1(1:4) fitParameters_picco2(1:3)]';
        cv_est_parGV = [cv_est_parGV_picco1 cv_est_parGV_picco2];
    else
        % Just use first pass parameters if no recirculation
        fitParameters = fitParameters_picco1(1:4);
        cv_est_parGV = cv_est_parGV_picco1;
    end
    
    % Store fit results in output structure
    AIF.fit.pesi = pesi;           % Fitting weights used
    AIF.fit.parameters = fitParameters;  % Fitted parameters
    AIF.fit.cv_est_parGV = cv_est_parGV;  % Parameter coefficient of variation
    
    % Calculate final fitted curve
    if options.aif.ricircolo
        AIF.fit.gv = GVfunction(fitParameters,options);  % With recirculation
    else
        AIF.fit.gv = GVfunction_picco1(fitParameters,options);  % First pass only
    end
    
    AIF.fit.time = options.time;  % Time points
    AIF.conc = AIFconc;          % Original concentration data
    
    
    % Store data for later plotting
    aif_slice_plot_data = struct();
    aif_slice_plot_data.image = image.img;
    aif_slice_plot_data.bound = image.bound;
    aif_slice_plot_data.ROI_x = AIF.ROI.x;
    aif_slice_plot_data.ROI_y = AIF.ROI.y;
    aif_slice_plot_data.voxels = AIF.voxels;
    aif_slice_plot_data.time = options.time;
    aif_slice_plot_data.conc = AIF.conc;
    aif_slice_plot_data.fit_gv = AIF.fit.gv;
    
    % Clean up figures if in verbose mode
    if options.display > 2
        handles = fieldnames(hf);
        for i=1:size(handles,1)
            try
                eval(['set(hf.' handles{i,:} ', ''CloseRequestFcn'', @my_closefcn);']);
            end
        end
    end
    
    end % extractAIF
    

    %% ------------------------------------------------------------------------
    function my_closefcn(src, ~)
        delete(src);
    end
    
    %%
    function optionsCopy = deepcopy(options)
        % Create a deep copy of the options structure
        optionsCopy = options;
        fields = fieldnames(options);
        for i = 1:numel(fields)
            if isstruct(options.(fields{i}))
                optionsCopy.(fields{i}) = deepcopy(options.(fields{i}));
            end
        end
    end
    
end % DSC_mri_aif (main function)