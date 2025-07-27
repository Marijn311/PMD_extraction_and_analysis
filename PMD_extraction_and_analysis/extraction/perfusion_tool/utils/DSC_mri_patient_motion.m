function [volumes]=DSC_mri_patient_motion(volumes, mask, options)
    
if options.display > 0
    disp('Detecting and correcting patient motion...')
end

% % if we detect to much patient motion, we will write the name of the file to the error_logs.txt file
% if motion > threshold
%     % if "error_logs.txt" does not exist, create it
%     if ~exist('error_logs.txt', 'file')
%         fileID = fopen('error_logs.txt', 'w');
%         fprintf(fileID, 'Error logs: \n');
%         fclose(fileID);
%     end
    
%     % write the name of the file to the error_logs.txt file
%     fileID = fopen('error_logs.txt', 'a');
%     fprintf(fileID, 'Too much motion detected for patient %s\n', options.perf_path);
%     fclose(fileID);
% end


    %Here you should implement a way to detect and or correct patient motion
    volumes = volumes; %return placeholder