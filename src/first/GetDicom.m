function [I, dicom_info] = GetDicom(file_name)
% Returns image and dicom_info for a certain file name. Also applies
% rescale parameters to the images. These parameters are available in the
% Phillips format, but not in the parsing made by OsiriX, so be careful.

I = [];
dicom_info = [];

try
    I = dicomread(file_name);
    dicom_info = dicominfo(file_name);
    I = double(I).*dicom_info.RescaleSlope + dicom_info.RescaleIntercept;
catch
    warning('File not found or other errors at reading DICOM file');
end


