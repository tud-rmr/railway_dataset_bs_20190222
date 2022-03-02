% Import reference data
%
%   Other m-files required: 
%       - setDatasetPaths
%       - loadFilenamePatterns
%       - importFromFile
%       - exportToFile
%   MAT-files required: none
%
%   See also: loadFilenamePatterns, processRawData

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 05-Jun-2020

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

% Load scripts ____________________________________________________________

loadFilenamePatterns

%% Checks

if ~exist('import_reference_data_executed','var')
    import_reference_data_executed = false;
elseif import_reference_data_executed
    fprintf('Import reference data script already executed!\n');
    return  
end % if

%% Processed Data
%
%   Filepath or filename has to contain session information, i.e. it has to 
%   contain a string which can be found with the regular 
%   expression '[Ss]ession\d{2}'. Valid strings are, 
%   e.g. 'Session01' or 'session01'.
%

ref_inatm200stn_data_paths = { ... 
                               fullfile('ref_inatm200stn_internal_ekf_data_session01.csv') ... 
%                                fullfile('2019-02-22_Session01','03_reference','ref_inatm200stn_session01.csv') ... 
                             };

ref_ekfFusion_javad_xsens_data_paths = { ... 
                                         fullfile('ref_ekfFusion_javad_xsens_data_session01.csv') ... 
%                                          fullfile('2019-02-22_Session01','03_reference','ref_ekfFusion_javad_xsens_session01.csv') ...
                                     };
                              
ref_ekfFusion_inatm200stn_data_paths = { ... 
                                         fullfile('ref_ekfFusion_inatm200stn_data_session01.csv') ... 
%                                          fullfile('2019-02-22_Session01','03_reference','ref_ekfFusion_inatm200stn_session01.csv') ...
                                       };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference: iNat M200 STN (internal EKF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(ref_inatm200stn_data_paths,ref_inatm200stn_reference_data_root_string,'DataType','reference');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[ref_inatm200stn_reference_data_root_string,'_session',sprintf('%02i',i)],ref_inatm200stn_reference_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',ref_inatm200stn_reference_data_root_string,out_data);
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference: EKF Fusion (Javad and Xsens)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(ref_ekfFusion_javad_xsens_data_paths,ref_ekfFusion_javad_xsens_reference_data_root_string,'DataType','reference');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[ref_ekfFusion_javad_xsens_reference_data_root_string,'_session',sprintf('%02i',i)],ref_ekfFusion_javad_xsens_reference_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',ref_ekfFusion_javad_xsens_reference_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference: EKF Fusion (iNat M200 STN GNSS and IMU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(ref_ekfFusion_inatm200stn_data_paths,ref_ekfFusion_inatm200stn_reference_data_root_string,'DataType','reference');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[ref_ekfFusion_inatm200stn_reference_data_root_string,'_session',sprintf('%02i',i)],ref_ekfFusion_inatm200stn_reference_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',ref_ekfFusion_inatm200stn_reference_data_root_string,out_data); 
clear out_data

%% Finish script

import_reference_data_executed = true;
