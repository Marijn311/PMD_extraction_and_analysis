function [mtt]=DSC_mri_mtt(cbv,cbf,options)
    
    if options.display > 0
        disp('   MTT');
    end

    % MTT is calculated as the ratio of CBV to CBF
    deconv_method= fieldnames(cbf);
    nr_deconv_methods = size(deconv_method, 1);


    for i = 1:nr_deconv_methods

        method = deconv_method{i};
        
        % Get the cbv and cbf maps
        cbf_map = cbf.(method).map;
        cbv_map = cbv;

        % Calculate the MTT map
        mtt_map = cbv_map ./ cbf_map;

        % Everything outside the mask seems to be set to nan, since outside the mask the cbf is zero and when we divide by zero we get nan
        % Convert any NaN values to 0
        mtt_map(isnan(mtt_map)) = 0;

        %assign to correct output
        mtt.(method) = mtt_map;


   

end