
% TTP is the time it takes for the signal intensity in a tissue voxel to reach its lowest point 
% (Which is the same as the time it takes to reach peak contrast agent concentration) during the DSC bolus passage.

% This metric depends on how the start of the recording relates to the start of contradst agent injection.
% Has this been standardized in some way? Does the recordig start at injection time, or a fixed x seconds before injection?
% If not than it makes no sense to compare TTP between different recordings and the 2 statistically significant regions i found might have been flukes

% If I want to fix it i need to determine either 
% A) confirm there is a fixed relative between start of recording and injection time. in this case no action required
% 
% B) I need define the start of the bolus myself. This can not be done for each pixel, the signal is too noisy for that.
%    I could define the start of the bolus when I average it over the entire brain.
%    This would mean that the ttp would still be higher for regions that are further away from the injection site, 
%    but at least it would be more comparable between recordings. 
%    There could still be relative difference between patient groups in the same recording that indicate a stroke. 
% 
% C) I could also make a rTTP by dividing the TTP by the TTP in a reference region, similar to how it is done with CBF and CBV.

% I will go with option C for now, because that is the same way the other metrics are calculated.

function [ttp]=DSC_mri_ttp(conc,mask,atlas,options)
    [~,ttp]=max(conc,[],4);
    ttp(not(logical(mask))) = 0;
   
end

