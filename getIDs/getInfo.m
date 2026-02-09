function [info,info_fullname,data_path] = getInfo(dataSet)

switch dataSet
    case 'Camk2'
        data_path = 'C:\Users\brndn\Downloads\CaMK2-ChR2\';
        info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\CaMK2-ChR2\';
        info_name = 'camk2_chr2_info_07102024.mat';
    case 'Thy1'
        data_path = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\';
        info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\';
        info_name = 'thy1_chr2_info_032822.mat';
    case 'PV Transgenic'
        data_path = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\';
        info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\PV-ChR2 Transgenic\';
        info_name = 'pv_chr2_transgenic_info_030323.mat';
    case 'PV Viral'
        data_path = 'C:\Users\brndn\Downloads\PV-ChR2\';
        info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\pv_chr2\';
        info_name = 'pv_chr2_info_060122.mat';
end
info_fullname = fullfile(info_path,info_name);
load(info_fullname,'info') % load info struct

end