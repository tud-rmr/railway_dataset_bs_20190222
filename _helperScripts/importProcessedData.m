% Import processed data
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

if ~exist('import_processed_data_executed','var')
    import_processed_data_executed = false;
elseif import_processed_data_executed
    fprintf('Import processed data script already executed!\n');
    return  
end % if

%% Processed Data
%
%   Filepath or filename has to contain session information, i.e. it has to 
%   contain a string which can be found with the regular 
%   expression '[Ss]ession\d{2}'. Valid strings are, 
%   e.g. 'Session01','session01','session02', etc
%

gnss_javad_processed_data_paths = { ... 
                                        fullfile('gnss_javad_processed_data_session01.csv') ... 
%                                         fullfile('2019-02-22_Session01','02_processed','gnss_javad_processed_data_session01.csv') ... 
                                      };
                      
imu_xsens_processed_data_paths = { ... 
                                       fullfile('imu_xsens_processed_data_session01.csv') ... 
%                                        fullfile('2019-02-22_Session01','02_processed','imu_xsens_processed_data_session01.csv') ... 
                                     };

speedDist_odometer_processed_data_paths = { ... 
                                            fullfile('speedDist_odometer_processed_data_session01.csv') ... 
%                                             fullfile('2019-02-22_Session01','02_processed','speedDist_odometer_processed_data_session01.csv') ... 
                                         };

speedDist_correvit_processed_data_paths = { ... 
                                                  fullfile('speedDist_correvit_processed_data_session01.csv') ... 
%                                                   fullfile('2019-02-22_Session01','02_processed','speedDist_correvit_processed_data_session01.csv') ... 
                                                };
                                 
speedDist_siemens_processed_data_paths = { ... 
                                               fullfile('speedDist_siemens_processed_data_session01.csv') ... 
%                                                fullfile('2019-02-22_Session01','02_processed','speedDist_siemens_processed_data_session01.csv') ... 
                                             };
                              
gnss_inatm200stn_processed_data_paths = { ... 
                                          fullfile('gnss_inatm200stn_processed_data_session01.csv') ... 
%                                           fullfile('2019-02-22_Session01','02_processed','gnss_inatm200stn_processed_data_session01.csv') ... 
                                        };
                         
imu_inatm200stn_processed_data_paths = { ... 
                                         fullfile('imu_inatm200stn_processed_data_session01.csv') ... 
%                                          fullfile('2019-02-22_Session01','02_processed','imu_inatm200stn_processed_data_session01.csv') ... 
                                       };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GNSS: Javad 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(gnss_javad_processed_data_paths,gnss_javad_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_javad_processed_data_root_string,'_session',sprintf('%02i',i)],gnss_javad_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_javad_processed_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU: Xsens 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(imu_xsens_processed_data_paths,imu_xsens_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[imu_xsens_processed_data_root_string,'_session',sprintf('%02i',i)],imu_xsens_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',imu_xsens_processed_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed/Distance: Odometer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(speedDist_odometer_processed_data_paths,speedDist_odometer_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[speedDist_odometer_processed_data_root_string,'_session',sprintf('%02i',i)],speedDist_odometer_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',speedDist_odometer_processed_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed/Distance: Correvit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(speedDist_correvit_processed_data_paths,speedDist_correvit_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[speedDist_correvit_processed_data_root_string,'_session',sprintf('%02i',i)],speedDist_correvit_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',speedDist_correvit_processed_data_root_string,out_data);
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed/Distance: Siemens Radar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(speedDist_siemens_processed_data_paths,speedDist_siemens_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[speedDist_siemens_processed_data_root_string,'_session',sprintf('%02i',i)],speedDist_siemens_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',speedDist_siemens_processed_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNSS: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(gnss_inatm200stn_processed_data_paths,gnss_inatm200stn_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_inatm200stn_processed_data_root_string,'_session',sprintf('%02i',i)],gnss_inatm200stn_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_inatm200stn_processed_data_root_string,out_data);
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(imu_inatm200stn_processed_data_paths,imu_inatm200stn_processed_data_root_string,'DataType','processed');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[imu_inatm200stn_processed_data_root_string,'_session',sprintf('%02i',i)],imu_inatm200stn_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',imu_inatm200stn_processed_data_root_string,out_data); 
clear out_data

%% Finish script

import_processed_data_executed = true;
