% Actin_kibra_xcorr.m
% Paper: Hino, Kapoor, Gubbala, Hannezo & Heisenberg
%        "Apical domain mechanosensation regulates tissue tension homeostasis"
%        Institute of Science and Technology Austria (ISTA)
%
% Computes spatial cross-correlation between Actin and Kibra fluorescence
% signals within segmented EVL cells. See image_analysis/README.md for details.

%% input
clear

Actin_img_name = 'Actin.tif';
Kibra_img_name = 'Kibra.tif';
Cellpose_img_name = 'Reference_surface_max_off2_delta5_th5_cp_masks.png';
Output_file_name = 'spatial_xcorr_results.xlsx';


Pixel_size = 0.1560; %um/pixel
Min_cell_area = 10; %pixel number
Erode_iteration = 6;

%%

Actin_img = double(imread(Actin_img_name));
Kibra_img = double(imread(Kibra_img_name));
Cellpose_img = double(imread(Cellpose_img_name));

Cell_ID = unique(Cellpose_img(:));
Cell_ID(Cell_ID == 0) = [];

Cell_number = numel(Cell_ID);

xcorrMap = cell(1,Cell_number);
Radius_um_binned = cell(1,Cell_number);
Radial_mean_xcorr = cell(1,Cell_number);
Nearest_peak_radius_after_mean = cell(1, Cell_number);
V_names_1 = cell(1,Cell_number+1);
V_names_2 = cell(1,Cell_number);

for i = 1:Cell_number
    target_cell = Cell_ID(i);
    Mask = (Cellpose_img == target_cell);
    Area_pixel = nnz(Mask);
    if Area_pixel < Min_cell_area
        xcorrMaps{i} = [];
        Radius_um_binned{i} = [];
        Radial_mean_xcorr{i} = [];
        Nearest_peak_radius_after_mean{i} = [];
        continue
    end
    
    Eroded_mask = Mask;
    se = strel('disk', 1, 0);
    for j = 1:Erode_iteration
        Eroded_mask = imerode(Eroded_mask, se);
    end
    Eroded_area_pixel = nnz(Eroded_mask);
    if Eroded_area_pixel < Min_cell_area
        xcorrMaps{i} = [];
        Radius_um_binned{i} = [];
        Radial_mean_xcorr{i} = [];
        Nearest_peak_radius_after_mean{i} = [];
        continue
    end
    [rr, cc] = find(Eroded_mask);
    rMin = min(rr);
    rMax = max(rr);
    cMin = min(cc);
    cMax = max(cc);

            cropMask = Eroded_mask(rMin:rMax, cMin:cMax);
            crop_Act = Actin_img(rMin:rMax, cMin:cMax);
            crop_Kib = Kibra_img(rMin:rMax, cMin:cMax);

            Act_pix = crop_Act(cropMask);
            Kib_pix = crop_Kib(cropMask);

            meanPix1 = mean(Act_pix);
            meanPix2 = mean(Kib_pix);

            Norm_Act = zeros(size(crop_Act));
            Norm_Kib = zeros(size(crop_Kib));
            Norm_Act(cropMask) = Act_pix - meanPix1;
            Norm_Kib(cropMask) = Kib_pix - meanPix2;
            
            denom = sqrt(sum(Norm_Act(:).^2) * sum(Norm_Kib(:).^2));
            corrMap = xcorr2(Norm_Act, Norm_Kib) ./ denom;
            xcorrMaps{i} = corrMap;
 
            
 %% Detect shift - nearest maximum to 0
            Center_Row = size(Norm_Act,1);
            Center_Col = size(Norm_Act,2);
            [Y_map, X_map] = ndgrid(1:size(corrMap,1), 1:size(corrMap,2));
            Y_shift_map = Y_map - Center_Row;
            X_shift_map = X_map - Center_Col;
            Y_shift_map_um = Y_shift_map * Pixel_size;
            X_shift_map_um = X_shift_map * Pixel_size;
            R_map_um = sqrt(Y_shift_map_um.^2 + X_shift_map_um.^2);
            
            max_R_um = min(size(Norm_Act,1),size(Norm_Act,2)) * Pixel_size;
            nBins = floor(max_R_um / Pixel_size);
            Radius_um = (0:nBins)' * Pixel_size;
            radial_mean = NaN(numel(Radius_um), 1);
            binIndex_map = round(R_map_um / Pixel_size);
            for b = 0:nBins
                radial_mean(b+1) = mean(corrMap(binIndex_map == b));
            end
           
 %% Detect shit - nearest local maximum to 0
            [peakValue, peakIdx] = findpeaks(radial_mean, Radius_um);
            peakIdx = min(peakIdx); 
            
        Radius_um_binned{i} = Radius_um;
        Radial_mean_xcorr{i} = radial_mean;
        Nearest_peak_radius_after_mean{i} = peakIdx;

end

lengths = cellfun(@numel, Radius_um_binned);
[Max_length, idx] = max(lengths);

Output_xcorr_data = NaN(Max_length, 1+Cell_number);
Output_after_mean_local_peak_data = NaN(1, Cell_number);

Output_xcorr_data(1:Max_length,1) = Radius_um_binned{idx};
V_names_1{1} = 'Distance_um';
for i = 1:Cell_number
    target_cell = Cell_ID(i);
    length_each_cell = numel(Radial_mean_xcorr{i});
    Output_xcorr_data(1:length_each_cell, i+1) = Radial_mean_xcorr{i};
    Output_after_mean_local_peak_data(1, i) = Nearest_peak_radius_after_mean{i};
    V_names_1{i+1} = sprintf('Cell_%d',target_cell);
    V_names_2{i} = sprintf('Cell_%d',target_cell);
end

Output_xcorr = array2table(Output_xcorr_data, 'VariableNames', V_names_1);
Output_after_mean_local_peak = array2table(Output_after_mean_local_peak_data, 'VariableNames', V_names_2);

writetable(Output_xcorr, Output_file_name, 'Sheet', 'radial_mean_xcorr');
writetable(Output_after_mean_local_peak, Output_file_name, 'Sheet', 'local_peak_distance_after_mean');