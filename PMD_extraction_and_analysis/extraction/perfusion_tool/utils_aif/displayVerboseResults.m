function displayVerboseResults(options, totalCandidates, ROI, ROIreg, hf, AIFslice, image, xROI, centro, ROIiniziale, AUC, threshold_auc, TTP, threshold_ttp, REG, threshold_reg, nT, nC, nR);

% Display results if verbose mode enabled
if options.display > 2
    disp(' ')
    disp(' Candidate voxel selection via regularity criteria')
    disp(['  Initial voxel count: ' num2str(totalCandidates)])
    disp(['  Surviving voxels: ' num2str(sum(sum(ROI)))])
    
    % Plot spatial map of kept/rejected voxels
    figure(hf.aif);
    subplot(2,3,3)
    imagesc(ROIreg)
    title('Regularity criteria - survived voxels') 
    set(gca,'xtick',[],'ytick',[])
    axis square
    
    % Plot time curves for kept and rejected voxels
    subplot(236)
    plot(options.time,zeros(1,nT),'b-')
    hold on
    plot(options.time,zeros(1,nT),'r-')
    plot(options.time,zeros(1,nT),'k-')
    for c=1:nC
        for r=1:nR
            if ROIreg(r,c)==1
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'r-')
            elseif ROIreg(r,c)==2
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'b-')
            end
        end
    end
    set(gca,'FontSize',12)
    legend('Accepted','Rejected')
    title('Irregularity Index')
    xlabel('time')
end


