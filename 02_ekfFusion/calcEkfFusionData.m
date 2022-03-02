% IMU/GNSS sensor-fusion positioning with EKF
%
%   Other m-files required: 
%       - setDatasetPaths
%       - loadFilenamePatterns
%       - importProcessedData
%       - runEkfFusion
%       - calcDiscreteKf
%       - stateTransitionFcn
%       - measurementFcn
%   MAT-files required*:
%       - processed GNSS data 
%       - processed IMU data 
%       * Run 'processRawData.m' to generate suitable data from the folder 
%         '[date]_[session##]/02_processed'.
%
%   See also: processRawData, runEkfFusion

%   Author: Hanno Winter
%   Date: 18-Feb-2020; Last revision: 12-Nov-2020

%% Settings (edit)

Tend = []; % duration of EKF fusion (if empty all data will be used)
Ts = [0.1]; % sample time of EKF fusion (if empty natural rate of data will be used)
session_selector = []; % session selector (if empty all available sessions will be used)

%% Init

% Set path ________________________________________________________________

saved_pwd = pwd;
pwd_parts = strfind(pwd,filesep);

if ~exist('dataset_paths_set','var') || (dataset_paths_set == false)
    
    dataset_paths_set = false;
    
    for i = 1:length(pwd_parts)
        
        if exist('setDatasetPaths','file')
            setDatasetPaths();
            dataset_paths_set = true;
            break;
        elseif exist('_fcns','dir')
            cd('_fcns');
        else
            cd('..')
        end % if
        
    end % for i
    
    cd(saved_pwd);
    
end % if

if ~dataset_paths_set
    error('Couldn''t add dataset folders to path!')
end % if

% Prepare calculatons _____________________________________________________

loadFilenamePatterns
importProcessedData

%% EKF Fusion: Javad and Xsens data

[x_hat_javad_xsens_ts,P_javad_xsens_ts,~,~,utm_zone] = runEkfFusion( ... 
                                                                     eval(gnss_javad_processed_data_root_string), ... 
                                                                     eval(imu_xsens_processed_data_root_string), ... 
                                                                     session_selector,... 
                                                                     Ts, ... 
                                                                     Tend ... 
                                                                   );
ref_ekfFusion_javad_xsens_processed_data = processEkfFusionData( ... 
                                                                 ref_ekfFusion_javad_xsens_reference_data_root_string, ... 
                                                                 ref_ekfFusion_javad_xsens_reference_data_root_string, ... 
                                                                 x_hat_javad_xsens_ts, ... 
                                                                 P_javad_xsens_ts, ... 
                                                                 utm_zone ... 
                                                               );
assignin('base',ref_ekfFusion_javad_xsens_reference_data_root_string,ref_ekfFusion_javad_xsens_processed_data);

%% EKF Fusion: iNat M200 STN GNSS and IMU data

[x_hat_inatm200stn_ts,P_inatm200stn_ts,~,~,utm_zone] = runEkfFusion( ... 
                                                                     eval(gnss_inatm200stn_processed_data_root_string), ... 
                                                                     eval(imu_inatm200stn_processed_data_root_string), ... 
                                                                     session_selector,... 
                                                                     Ts, ... 
                                                                     Tend ... 
                                                                   );
ref_ekfFusion_inatm200stn_processed_data = processEkfFusionData( ... 
                                                                 ref_ekfFusion_inatm200stn_reference_data_root_string, ... 
                                                                 ref_ekfFusion_inatm200stn_reference_data_root_string, ... 
                                                                 x_hat_inatm200stn_ts, ... 
                                                                 P_inatm200stn_ts, ... 
                                                                 utm_zone ... 
                                                               );
assignin('base',ref_ekfFusion_inatm200stn_reference_data_root_string,ref_ekfFusion_inatm200stn_processed_data);
