% Plot processed output data of simulink EKF Fusion model
%
%   Other m-files required: 
%       * setDatasetPaths
%       * importReferenceData
%       * importProcessedData
%       * getLatLonErrorEllipsePoints
%   MAT-files required: none
%
%   See also: loadFilenamePatterns, importReferenceData, 
%             importProcessedData, processEkfFusionData

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 24-Nov-2020

%% Init

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

%% EKF positioning

if(1)
    session_selector = [1];
    
    plot_error_ellipses = 0; % enable or disable error-ellipse plots
    ts_error_ellipses_javed_xsens = 1; % distance between error-ellipse plots in sec
    ts_error_ellipses_inatm200stn = 1; % distance between error-ellipse plots in sec
    max_major_error_ellipses = 10; % maximum major semi-axis size of error ellipses to plot
    
    % Pre-Calculations ____________________________________________________
    importReferenceData
    importMapData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        
        if plot_error_ellipses
            error_ellipses_inatm200stn_indices = 1:ceil(ts_error_ellipses_inatm200stn/seconds(mean(diff(ref_ekfFusion_inatm200stn_data{session_i}.Time)))):length(ref_ekfFusion_inatm200stn_data{session_i}.Time);
            error_ellipses_inatm200stn_indices = error_ellipses_inatm200stn_indices(ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseMajor_m(error_ellipses_inatm200stn_indices)' < max_major_error_ellipses);
            [~,ref_ekfFusion_inatm200stn_error_ellipse] = ... 
                getLatLonErrorEllipsePoints( ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}{error_ellipses_inatm200stn_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseMajor_m(error_ellipses_inatm200stn_indices)*2.45, ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseMinor_m(error_ellipses_inatm200stn_indices)*2.45, ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseOrientation_deg(error_ellipses_inatm200stn_indices), ... 
                                             360 ... 
                                           );
            error_ellipses_javad_xsens_indices = 1:ceil(ts_error_ellipses_javed_xsens/seconds(mean(diff(ref_ekfFusion_javad_xsens_data{session_i}.Time)))):length(ref_ekfFusion_javad_xsens_data{session_i}.Time);
            error_ellipses_javad_xsens_indices = error_ellipses_javad_xsens_indices(ref_ekfFusion_javad_xsens_data{session_i}.ErrorEllipseMajor_m(error_ellipses_javad_xsens_indices)' < max_major_error_ellipses);
            [~,ref_ekfFusion_javad_xsens_error_ellipse] = ... 
                getLatLonErrorEllipsePoints( ... 
                                             ref_ekfFusion_javad_xsens_data{session_i}{error_ellipses_javad_xsens_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                             ref_ekfFusion_javad_xsens_data{session_i}.ErrorEllipseMajor_m(error_ellipses_javad_xsens_indices)*2.45, ... 
                                             ref_ekfFusion_javad_xsens_data{session_i}.ErrorEllipseMinor_m(error_ellipses_javad_xsens_indices)*2.45, ... 
                                             ref_ekfFusion_javad_xsens_data{session_i}.ErrorEllipseOrientation_deg(error_ellipses_javad_xsens_indices), ... 
                                             360 ... 
                                           );
        end % if
        
        % Plot ____________________________________________________________
        figure_name = ['EKF Postions (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name);
        clear h_plot    
        h_plot = gobjects(0);
        
        if verLessThan('matlab', '9.5') 
            
            hold on; grid on;
            [h_tracks, h_map] = plotXyTrackMap(map_xy_track_map); 
            set(h_tracks(:),'Visible','off');
            set(h_map,'LineWidth',2,'Color','k','LineStyle','-','DisplayName','Track Map');

            h_plot(end+1) = plot(ref_inatm200stn_internal_ekf_data{session_i}.Longitude_deg,ref_inatm200stn_internal_ekf_data{session_i}.Latitude_deg,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
            h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_data{session_i}.Longitude_deg,ref_ekfFusion_inatm200stn_data{session_i}.Latitude_deg,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
            h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_data{session_i}.Longitude_deg,ref_ekfFusion_javad_xsens_data{session_i}.Latitude_deg,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');
            if plot_error_ellipses
                h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_error_ellipse(2,:),ref_ekfFusion_inatm200stn_error_ellipse(1,:),'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF 2.45\sigma-Error-Ellipse (\approx95%) (iNatM200Stn GNSS+IMU)');
                h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_error_ellipse(2,:),ref_ekfFusion_javad_xsens_error_ellipse(1,:),'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF 2.45\sigma-Error-Ellipse (\approx95%) (Javad+Xsens)');
            end % if
            
            % You can try
            % https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
            % to create satellite plots with Google Maps API
            if exist('plot_google_map')            
                plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
            end % if
            
            h_legend = legend([h_map,h_plot]);
            set(h_legend,'Location','best')
            xlabel('longitude [deg]')
            ylabel('latitude [deg]')
            
        else
            
            [h_tracks, h_map] = plotXyTrackMap(map_xy_track_map); 
            set(h_tracks(:),'Visible','off');
            set(h_map,'LineWidth',2,'Color','k','LineStyle','-','DisplayName','Track Map');

            h_plot(end+1) = geoplot(ref_inatm200stn_internal_ekf_data{session_i}.Latitude_deg,ref_inatm200stn_internal_ekf_data{session_i}.Longitude_deg,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
            h_plot(end+1) = geoplot(ref_ekfFusion_inatm200stn_data{session_i}.Latitude_deg,ref_ekfFusion_inatm200stn_data{session_i}.Longitude_deg,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)'); hold on;   
            h_plot(end+1) = geoplot(ref_ekfFusion_javad_xsens_data{session_i}.Latitude_deg,ref_ekfFusion_javad_xsens_data{session_i}.Longitude_deg,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');
            if plot_error_ellipses
                h_plot(end+1) = geoplot(ref_ekfFusion_inatm200stn_error_ellipse(1,:),ref_ekfFusion_inatm200stn_error_ellipse(2,:),'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF 2.45\sigma-Error-Ellipse (\approx95%) (iNatM200Stn GNSS+IMU)');
                h_plot(end+1) = geoplot(ref_ekfFusion_javad_xsens_error_ellipse(1,:),ref_ekfFusion_javad_xsens_error_ellipse(2,:),'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF 2.45\sigma-Error-Ellipse (\approx95%) (Javad+Xsens)');
            end % if
            h_legend = legend([h_map,h_plot]);
            set(h_legend,'Location','best')
            geobasemap('topographic'); % 'satellite','topographic','none'
           
        end        
        
    end % for session_i
end % if

%% EKF motion data

if(1)
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importReferenceData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        ref_inatm200stn_internal_ekf_utc_time = datetime(seconds(ref_inatm200stn_internal_ekf_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_ekfFusion_inatm200stn_utc_time = datetime(seconds(ref_ekfFusion_inatm200stn_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_ekfFusion_javad_xsens_utc_time = datetime(seconds(ref_ekfFusion_javad_xsens_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['EKF Motion data (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.DistanceVehicle_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
        h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_utc_time,ref_ekfFusion_javad_xsens_data{session_i}.DistanceVehicle_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('d_{vehicle} [m]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.VelocityVehicle_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
        h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_utc_time,ref_ekfFusion_javad_xsens_data{session_i}.VelocityVehicle_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('v_{vehicle} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(ref_inatm200stn_internal_ekf_utc_time,ref_inatm200stn_internal_ekf_data{session_i}.AccX_mss,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
        h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.AccX_mss,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
        h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_utc_time,ref_ekfFusion_javad_xsens_data{session_i}.AccX_mss,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');   
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('a_{vehicle} [m/s^2]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if

%% Distance and Speed Data

if(1)  
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    importReferenceData        
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        gnss_inatm200stn_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_javad_utc_time = datetime(seconds(gnss_javad_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        speedDist_odometer_utc_time = datetime(seconds(speedDist_odometer_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        speedDist_correvit_utc_time = datetime(seconds(speedDist_correvit_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        speedDist_siemens_utc_time = datetime(seconds(speedDist_siemens_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_inatm200stn_internal_ekf_utc_time = datetime(seconds(ref_inatm200stn_internal_ekf_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_ekfFusion_inatm200stn_utc_time = datetime(seconds(ref_ekfFusion_inatm200stn_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_ekfFusion_javad_xsens_utc_time = datetime(seconds(ref_ekfFusion_javad_xsens_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['EKF Vehicle Distances and Speeds (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.DistanceVehicle_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
        h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_utc_time,ref_ekfFusion_javad_xsens_data{session_i}.DistanceVehicle_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.DistanceVehicle_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.DistanceVehicle_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.DistanceVehicle_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('d_{vehicle} [m]')

        clear h_plot    
        h_plot = gobjects(0);   
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_inatm200stn_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityGround_ms,'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn GNSS');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.VelocityGround_ms,'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad GNSS');
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_plot(end+1) = plot(ref_inatm200stn_internal_ekf_utc_time,ref_inatm200stn_internal_ekf_data{session_i}.VelocityGround_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
        h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.VelocityGround_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
        h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_utc_time,ref_ekfFusion_javad_xsens_data{session_i}.VelocityGround_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');     
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('v_{ground} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_inatm200stn_utc_time,gnss_inatm200stn_processed_data{session_i}.Heading_deg,'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn GNSS');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.Heading_deg,'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad GNSS');
        h_plot(end+1) = plot(ref_inatm200stn_internal_ekf_utc_time,ref_inatm200stn_internal_ekf_data{session_i}.Heading_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
        h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.Heading_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (iNatM200Stn GNSS+IMU)');
        h_plot(end+1) = plot(ref_ekfFusion_javad_xsens_utc_time,ref_ekfFusion_javad_xsens_data{session_i}.Heading_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF (Javad+Xsens)');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('heading angle [deg]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if

%% Helper Functions

function [h_tracks, h_map] = plotXyTrackMap(xy_tack_map)
% [h_tracks, h_map] = plotXyTrackMap(xy_tack_map)
% 
%   Create plot of (X,Y)-track-map
%
%   In:
%       xy_track_map     (X,Y)-track-map data as table
%
%   Out:    
%       h_tracks         Array of line objects for each track
%       h_map            Line objects including the whole map
% 
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: none

%   Author: Hanno Winter
%   Date: 24-Mar-2020; Last revision: 13-May-2020

%% Init

tracks = unique(xy_tack_map(:,1).track_num);
h_tracks = gobjects(length(tracks),1);
map_data = table();

%% Calculations

for i = 1:length(tracks)
    
    track_i = tracks(i);    
    track_i_selector = (xy_tack_map(:,1).track_num == track_i);    
    track_i_data = xy_tack_map(track_i_selector,:);
    
    track_i_map_data = track_i_data;
    track_i_map_data{end+1,:} = nan;
    map_data = [map_data; track_i_map_data];
    
    if verLessThan('matlab', '9.5')     
        h_tracks(i) = plot(track_i_data.lon,track_i_data.lat,'DisplayName',sprintf('Track: %i',track_i)); hold on;
    else
        h_tracks(i) = geoplot(track_i_data.lat,track_i_data.lon,'DisplayName',sprintf('Track: %i',track_i)); hold on;
        geobasemap('topographic')        
    end % if
    
end % for i

% Plot whole map at once
if verLessThan('matlab', '9.5') 
    h_map = plot(map_data.lon,map_data.lat);
else
    h_map = geoplot(map_data.lat,map_data.lon);
end % if

end % function
