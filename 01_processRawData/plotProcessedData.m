% Plot processed raw data
%
%   Other m-files required: setDatasetPaths, importProcessedData
%   MAT-files required: none
%
%   See also: loadFilenamePatterns, importProcessedData

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

%% GNSS Positions

if(1)        
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________    
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Lat-Lon-Positions (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); 

        clear h_plot    
        h_plot = gobjects(0);  
        
        if verLessThan('matlab', '9.5') 
            % You can try
            % https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
            % to create satellite plots with Google Maps API
            hold on; grid on;
            h_plot(end+1) = plot(gnss_inatm200stn_processed_data{session_i}.Longitude_deg,gnss_inatm200stn_processed_data{session_i}.Latitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
            h_plot(end+1) = plot(gnss_javad_processed_data{session_i}.Longitude_deg,gnss_javad_processed_data{session_i}.Latitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');    
            legend(h_plot);
            xlabel('longitude [deg]')
            ylabel('latitude [deg]')
        else
           h_plot(end+1) = geoplot(gnss_inatm200stn_processed_data{session_i}.Latitude_deg,gnss_inatm200stn_processed_data{session_i}.Longitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn'); hold on;
           h_plot(end+1) = geoplot(gnss_javad_processed_data{session_i}.Latitude_deg,gnss_javad_processed_data{session_i}.Longitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');     
           geobasemap('topographic');
           legend(h_plot);
        end
    
    end % for session_i
end % if

%% GNSS Altitudes

if(1)
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________    
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_javad_utc_time = datetime(seconds(gnss_javad_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Altitudes from Ellipsoid (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold on; grid on;

        clear h_plot    
        h_plot = gobjects(0);         
        
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.AltitudeEllipsoid_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.AltitudeEllipsoid_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');    

        legend(h_plot);
        xlabel('UTC Time')
        ylabel('height_{ellip.} [m]')
    
    end % for session_i
end % if

%% GNSS Speeds

if(1)
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_javad_utc_time = datetime(seconds(gnss_javad_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS NED-Velocities (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;    
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityNorth_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.VelocityNorth_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('v_{north} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityEast_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.VelocityEast_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('v_{east} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityDown_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.VelocityDown_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('v_{down} [m/s]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if

%% Ground Speeds (GNSS)

if(1)   
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_javad_utc_time = datetime(seconds(gnss_javad_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Speed over Ground (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(2,1,1); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('v_{ground} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(2,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.Heading_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.Heading_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('heading [deg]')

        linkaxes([ax1,ax2],'x');
        
    end % for session_i
end % if

%% Ground Distances and Speeds

if(1)    
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________    
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_javad_utc_time = datetime(seconds(gnss_javad_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        speedDist_odometer_utc_time = datetime(seconds(speedDist_odometer_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        speedDist_correvit_utc_time = datetime(seconds(speedDist_correvit_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        speedDist_siemens_utc_time = datetime(seconds(speedDist_siemens_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['Vehicle Distances and Speeds (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;   

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(5,1,1); hold on; grid on;
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.DistanceVehicle_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.DistanceVehicle_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.DistanceVehicle_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('d_{vehicle} [m]')

        clear h_plot    
        h_plot = gobjects(0);   
        ax2 = subplot(5,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('v_{ground} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(5,1,3); hold on; grid on;
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.DirectionVehicle,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.DirectionVehicle,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.DirectionVehicle,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_legend = legend(h_plot);
        set(h_legend,'Location','northeast')
        xlabel('UTC Time')
        ylim([-2 2])
        yticks([-1 0 1])
        yticklabels({'backwards','unkown','forwards'})
        ylabel('direction')

        clear h_plot
        h_plot = gobjects(0); 
        ax4 = subplot(5,1,4); hold on; grid on;        
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.StandstillFlag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.StandstillFlag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.StandstillFlag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylim([-0.5 1.5])
        yticks([0 1])
        yticklabels({'unknown','standstill'})
        ylabel('standstill flag')
        
        clear h_plot
        h_plot = gobjects(0); 
        ax5 = subplot(5,1,5); hold on; grid on;
        h_plot(end+1) = plot(speedDist_odometer_utc_time,speedDist_odometer_processed_data{session_i}.MotionFlag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Odometer');
        h_plot(end+1) = plot(speedDist_correvit_utc_time,speedDist_correvit_processed_data{session_i}.MotionFlag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optical');
        h_plot(end+1) = plot(speedDist_siemens_utc_time,speedDist_siemens_processed_data{session_i}.MotionFlag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Siemens Radar');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylim([-0.5 1.5])
        yticks([0 1])
        yticklabels({'unknown','moving'})
        ylabel('motion flag')

        linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
        
    end % for session_i
end % if

%% GNSS Position Uncertainties

if(1)    
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_javad_utc_time = datetime(seconds(gnss_javad_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Position Uncertainties (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.LatitudeSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.LatitudeSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('\sigma_{latitue} [m]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.LongitudeSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.LongitudeSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('\sigma_{longitude} [m]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.AltitudeEllipsoidSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(gnss_javad_utc_time,gnss_javad_processed_data{session_i}.AltitudeEllipsoidSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Javad');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('\sigma_{height} [m]')

        linkaxes([ax1,ax2,ax3],'x');  
        
    end % for session_i
end % if

%% IMU Accelerations

if(1)    
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        imu_inatm200stn_utc_time = datetime(seconds(imu_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        imu_xsens_utc_time = datetime(seconds(imu_xsens_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['IMU Accelerations (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(imu_inatm200stn_utc_time,imu_inatm200stn_processed_data{session_i}.AccX_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(imu_xsens_utc_time,imu_xsens_processed_data{session_i}.AccX_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Xsens MTi-G-700'); 
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('a_x [m/s^2]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(imu_inatm200stn_utc_time,imu_inatm200stn_processed_data{session_i}.AccY_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(imu_xsens_utc_time,imu_xsens_processed_data{session_i}.AccY_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Xsens MTi-G-700');    
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('a_y [m/s^2]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(imu_inatm200stn_utc_time,imu_inatm200stn_processed_data{session_i}.AccZ_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(imu_xsens_utc_time,imu_xsens_processed_data{session_i}.AccZ_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Xsens MTi-G-700');    
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('a_z [m/s^2]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if

%% IMU Turn Rates

if(1)  
    session_selector = [1];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________        
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        imu_inatm200stn_utc_time = datetime(seconds(imu_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        imu_xsens_utc_time = datetime(seconds(imu_xsens_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['IMU Turn Rates (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(imu_inatm200stn_utc_time,imu_inatm200stn_processed_data{session_i}.TurnRateX_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(imu_xsens_utc_time,imu_xsens_processed_data{session_i}.TurnRateX_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Xsens MTi-G-700');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('w_x [deg/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(imu_inatm200stn_utc_time,imu_inatm200stn_processed_data{session_i}.TurnRateY_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(imu_xsens_utc_time,imu_xsens_processed_data{session_i}.TurnRateY_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Xsens MTi-G-700');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('w_y [deg/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(imu_inatm200stn_utc_time,imu_inatm200stn_processed_data{session_i}.TurnRateZ_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_plot(end+1) = plot(imu_xsens_utc_time,imu_xsens_processed_data{session_i}.TurnRateZ_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Xsens MTi-G-700');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC Time')
        ylabel('w_z [deg/s]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if
