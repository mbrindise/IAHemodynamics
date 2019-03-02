function [] = plot_vel_profiles(cadat,fslash,CREATE_MASKS,CREATE_FRES_DIS)
% Plots the velocity profiles for the cerebral aneurysm velocity fields
% The profiles plotted are at all inlets and outlets and any other place
% specified by the user. The inlet and outlet locations are plotted
% automatically based on the cross section data.
% CREATE_XCS - Logical that indicates whether the cross sections for
% plotting the velocity profiles needs to be specified by user or if they
% can be loaded. The cross sections for plotting are not necessarily the
% cross sections for computing flow rate - it depends on user inputs here
% CREATE_MASKS - Logical that indicates whether the masks for velocity
% distribution comparisons needs to be created.
% **Note - If CREATE_XCS and CREATE_MASKS are set to 0, it means the files
% needed for plotting have already been created and will be loaded.

%% LOAD DATA AND INPUTS
%%% IDENTIFY PEAK SYSTOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set file number of peak systole for each modality
% Peak systole occurs at t = 0.5265 seconds
% Peak diastole occurs at t = 0.324 seconds
% Mid diastole occurs at t = 0.081 seconds
% Mid systole occurs at t = 0.2277 seconds
if strcmp(cadat.FILE.anynum,'111')
    ps_mri = 14; % For ANY111 = 14, for ANY007 = 5
    ps_stb = 213; % For ANY111 = 211, for ANY007 = 102
    ps_cfd = 323; % For ANY111 = 323, for ANY007 = 102
    % Set the file numbers
    mri_file_s = 1; mri_file_e = 20; mri_file_skip = 1;
    stb_file_s = 1; stb_file_e = 304; stb_file_skip = 1;
    cfd_file_s = 1; cfd_file_e = 541; cfd_file_skip = 1;
else
    ps_mri = 4; % For ANY111 = 14, for ANY007 = 4
    ps_stb = 97; % For ANY111 = 211, for ANY007 = 102
    ps_cfd = 102; % For ANY111 = 323, for ANY007 = 102
    % Set the file numbers
    mri_file_s = 1; mri_file_e = 13; mri_file_skip = 1;
    stb_file_s = 1; stb_file_e = 365; stb_file_skip = 1;
    cfd_file_s = 1; cfd_file_e = 386; cfd_file_skip = 1;
end
    
%%% Load the velocity masks for each modality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the MRI velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_MRI_masked_velocity_grid.mat']);
mri_velmask = fp.velmask; % Velocity grid mask
mri_maskX = fp.xmask;
mri_maskY = fp.ymask;
mri_maskZ = fp.zmask;
% Load the STB velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_STB_masked_velocity_grid-voxel_averaged.mat']);
stb_velmask = fp.velmask; % Velocity grid mask
stb_maskX = fp.xmask;
stb_maskY = fp.ymask;
stb_maskZ = fp.zmask;
% Load the CFD velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_CFD_masked_velocity_grid-voxel_averaged.mat']);
cfd_velmask = fp.velmask; % Velocity grid mask
cfd_maskX = fp.xmask;
cfd_maskY = fp.ymask;
cfd_maskZ = fp.zmask;
% Load the full resolution STB velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_STB_masked_velocity_grid.mat']);
stb_fres_velmask = fp.velmask; % Velocity grid mask
stb_fres_maskX = fp.xmask;
stb_fres_maskY = fp.ymask;
stb_fres_maskZ = fp.zmask;
% Load the full resolution CFD velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_CFD_masked_velocity_grid.mat']);
cfd_fres_velmask = fp.velmask; % Velocity grid mask
cfd_fres_maskX = fp.xmask;
cfd_fres_maskY = fp.ymask;
cfd_fres_maskZ = fp.zmask;

%%% Load the STL masks for each modality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load peak systole velocity field for all modalities %%%%%%%%%%%%%%%%%%%
% Load MRI velocity field
if strcmp(cadat.FILE.anynum,'111')
    fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'MRI',fslash,'vel_registered_masked',fslash,'any',cadat.FILE.anynum,'_MRI_regmask_',num2str(ps_mri,'%05i'),'.mat']);
else
    fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'MRI',fslash,'vel_gridded',fslash,'any',cadat.FILE.anynum,'_MRI_regmask_gridded_',num2str(ps_mri,'%05i'),'.mat']);
end
mri_x = fp.x; mri_y = fp.y; mri_z = fp.z;
mri_u = fp.u; mri_v = fp.v; mri_w = fp.w;
% Load STB velocity field
fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'STB',fslash,'vel_p3_voxavg',fslash,'any',cadat.FILE.anynum,'_STB_v1_p3_voxavg_',num2str(ps_stb,'%05i'),'.mat']);
stb_x = fp.x; stb_y = fp.y; stb_z = fp.z;
stb_u = fp.u; stb_v = fp.v; stb_w = fp.w;
% Load CFD velocity field
fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'CFD',fslash,'vel_p3_voxavg',fslash,'any',cadat.FILE.anynum,'_CFD_p3_voxavg_',num2str(ps_cfd,'%05i'),'.mat']);
cfd_x = fp.x; cfd_y = fp.y; cfd_z = fp.z;
cfd_u = fp.u; cfd_v = fp.v; cfd_w = fp.w;

%%% Load flow rate data for all modalities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load MRI flow rate
fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'_MRI_flow_rate_data.mat']);
flowrate_mri = fp.flow_rate;
% Load STB flow rate
fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'_STB_flow_rate_data.mat']);
flowrate_stb = fp.flow_rate;
% Load CFD flow rate
fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'_CFD_flow_rate_data.mat']);
flowrate_cfd = fp.flow_rate;
% Load full resolution STB flow rate
fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'_STB_full_resolution_flow_rate_data.mat']);
flowrate_full_stb = fp.flow_rate;
% Load full resolution CFD flow rate
fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'_CFD_full_resolution_flow_rate_data.mat']);
flowrate_full_cfd = fp.flow_rate;

%%% Create mask for aneurysm velocity distribution comparison %%%%%%%%%%%%%
if CREATE_MASKS
    %%% MRI Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the MRI data and have the user change all points in mask to 1
    [mri_xi,mri_yi,mri_zi] = meshgrid(1:size(mri_velmask,2),1:size(mri_velmask,1),1:size(mri_velmask,3));
    % Reshape the points to line plots
    mri_xil = reshape(mri_xi,[size(mri_xi,1)*size(mri_xi,2)*size(mri_xi,3),1]);
    mri_yil = reshape(mri_yi,[size(mri_xi,1)*size(mri_xi,2)*size(mri_xi,3),1]);
    mri_zil = reshape(mri_zi,[size(mri_xi,1)*size(mri_xi,2)*size(mri_xi,3),1]);
    mri_vmaskl = reshape(mri_velmask,[size(mri_xi,1)*size(mri_xi,2)*size(mri_xi,3),1]);
    % limit the points to only those in the domain
    kp_inds = ~isnan(mri_vmaskl);
    mri_xik = mri_xil(kp_inds);
    mri_yik = mri_yil(kp_inds);
    mri_zik = mri_zil(kp_inds);
    figure(81); hold off; scatter3(mri_x,mri_y,mri_z,'k','MarkerFaceColor','k'); title('MRI'); view(2);
    xlabel('x'); ylabel('y'); zlabel('z');
    mri_any_mask = mri_velmask;
    fprintf('\nSet all points in the mri mask to 1 in the mri_any_mask array\n')
    keyboard
    % Ensure only points in the mask are used
    mri_any_mask_ref = mri_any_mask.*(~isnan(mri_velmask));
    % Reshape the points for scattering
    [Xm,Ym,Zm] = meshgrid(mri_maskX,mri_maskY,mri_maskZ);
    xmL = reshape(Xm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    ymL = reshape(Ym,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    zmL = reshape(Zm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    any_mask = reshape(mri_any_mask_ref,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    % Keep only the portions of the mask that equal 1
    kp_inds = (any_mask == 1);
    xmL = xmL(kp_inds==1);
    ymL = ymL(kp_inds==1);
    zmL = zmL(kp_inds==1);
    figure(81); hold on; scatter3(xmL,ymL,zmL,'r','MarkerFaceColor','r')
    fprintf('\nAdjust the mask points (xmL, ymL, zmL) if needed, after this step the mask will be saved\n')
    keyboard
    % Plot points after user adjusts Adjust points 
    figure(80); hold off; scatter3(mri_x,mri_y,mri_z,'k','MarkerFaceColor','k'); title('MRI'); view(2);
    hold on; scatter3(xmL,ymL,zmL,'r','MarkerFaceColor','r')
    % Save the MRI mask
    any_mask = [xmL,ymL,zmL];
    save([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_MRI_aneursym_mask.mat'],'any_mask')
    
    %%% STB Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the MRI data and have the user change all points in mask to 1
    figure(80); hold off; scatter3(stb_x,stb_y,stb_z,'k','MarkerFaceColor','k'); title('MRI'); view(2);
    stb_any_mask = stb_velmask;
    fprintf('\nSet all points in the stb mask to 1 in the stb_any_mask array\n')
    keyboard
    % Ensure only points in the mask are used
    stb_any_mask_ref = stb_any_mask.*(~isnan(stb_velmask));
    % Reshape the points for scattering
    [Xm,Ym,Zm] = meshgrid(stb_maskX,stb_maskY,stb_maskZ);
    xmL = reshape(Xm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    ymL = reshape(Ym,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    zmL = reshape(Zm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    any_mask = reshape(stb_any_mask_ref,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    % Keep only the portions of the mask that equal 1
    kp_inds = (any_mask == 1);
    xmL = xmL(kp_inds==1);
    ymL = ymL(kp_inds==1);
    zmL = zmL(kp_inds==1);
    figure(80); hold on; scatter3(xmL,ymL,zmL,'r','MarkerFaceColor','r')
    fprintf('\nAdjust the mask points (xmL, ymL, zmL) if needed, after this step the mask will be saved\n')
    keyboard
    % Plot points after user adjusts Adjust points 
    figure(80); hold off; scatter3(stb_x,stb_y,stb_z,'k','MarkerFaceColor','k'); title('STB'); view(2);
    hold on; scatter3(xmL,ymL,zmL,'r','MarkerFaceColor','r')
    % Save the STB mask
    any_mask = [xmL,ymL,zmL];
    save([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_STB_aneursym_mask.mat'],'any_mask')
    
    
    %%% CFD Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the CFD data and have the user change all points in mask to 1
    figure(80); hold off; scatter3(cfd_x,cfd_y,cfd_z,'k','MarkerFaceColor','k'); title('CFD'); view(2);
    cfd_any_mask = cfd_velmask;
    fprintf('\nSet all points in the cfd mask to 1 in the cfd_any_mask array\n')
    keyboard
    % Ensure only points in the mask are used
    cfd_any_mask_ref = cfd_any_mask.*(~isnan(cfd_velmask));
    % Reshape the points for scattering
    [Xm,Ym,Zm] = meshgrid(cfd_maskX,cfd_maskY,cfd_maskZ);
    xmL = reshape(Xm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    ymL = reshape(Ym,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    zmL = reshape(Zm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    any_mask = reshape(cfd_any_mask_ref,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
    % Keep only the portions of the mask that equal 1
    kp_inds = (any_mask == 1);
    xmL = xmL(kp_inds==1);
    ymL = ymL(kp_inds==1);
    zmL = zmL(kp_inds==1);
    figure(80); hold on; scatter3(xmL,ymL,zmL,'r','MarkerFaceColor','r')
    fprintf('\nAdjust the mask points (xmL, ymL, zmL) if needed, after this step the mask will be saved\n')
    keyboard
    % Plot points after user adjusts Adjust points 
    figure(80); hold off; scatter3(cfd_x,cfd_y,cfd_z,'k','MarkerFaceColor','k'); title('CFD'); view(2);
    hold on; scatter3(xmL,ymL,zmL,'r','MarkerFaceColor','r')
    % Save the CFD mask
    any_mask = [xmL,ymL,zmL];
    save([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_CFD_aneursym_mask.mat'],'any_mask')
else
    % Load the masks
    fp = load([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_MRI_aneursym_mask.mat']);
    mri_any_mask = fp.any_mask;
    fp = load([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_STB_aneursym_mask.mat']);
    stb_any_mask = fp.any_mask;
    fp = load([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_CFD_aneursym_mask.mat']);
    cfd_any_mask = fp.any_mask;
end

%%% Determine if any flow rates are negative, if so reverse them
q_stb = flowrate_stb.Q;
q_stb = q_stb*60000; % Convert to L/min
q_full_stb = flowrate_full_stb.Q;
q_full_stb = q_full_stb*60000; % Convert to L/min
q_mri = flowrate_mri.Q;
q_mri = q_mri*60000; % Convert to L/min
q_cfd = flowrate_cfd.Q;
q_cfd = q_cfd*60000; % Convert to L/min
q_full_cfd = flowrate_full_cfd.Q;
q_full_cfd = q_full_cfd*60000; % Convert to L/min
% Convert to mL/s
q_stb = q_stb*1000/60;
q_full_stb = q_full_stb*1000/60;
q_mri = q_mri*1000/60;
q_cfd = q_cfd*1000/60;
q_full_cfd = q_full_cfd*1000/60;

for zz = 1:1:size(q_stb,2)
    %q_stb(:,zz) = smooth(q_stb(:,zz),'moving',5);
    neg_perc = sum(q_stb(:,zz) < 0)/length(q_stb(:,zz));
    if neg_perc > 0.25, q_stb(:,zz) = -q_stb(:,zz); end
    neg_perc = sum(q_cfd(:,zz) < 0)/length(q_cfd(:,zz));
    if neg_perc > 0.25, q_cfd(:,zz) = -q_cfd(:,zz); end
    neg_perc = sum(q_full_stb(:,zz) < 0)/length(q_full_stb(:,zz));
    if neg_perc > 0.25, q_full_stb(:,zz) = -q_full_stb(:,zz); end
    neg_perc = sum(q_full_cfd(:,zz) < 0)/length(q_full_cfd(:,zz));
    if neg_perc > 0.25, q_full_cfd(:,zz) = -q_full_cfd(:,zz); end
    neg_perc = sum(q_mri(:,zz) < 0)/length(q_mri(:,zz));
    if neg_perc > 0.25, q_mri(:,zz) = -q_mri(:,zz); end
end

%%% Create the time vector for each modality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_mri = transpose(0:cadat.options.dt.mri:cadat.options.dt.mri*(size(q_mri,1)-1));
t_stb = transpose(0:cadat.options.dt.stb:cadat.options.dt.stb*(size(q_stb,1)-1));
t_cfd = transpose(0:cadat.options.dt.cfd:cadat.options.dt.cfd*(size(q_cfd,1)-1));
Tpc = max(t_mri); % Pulsatile cycle time
t_mri = t_mri/Tpc; % Non-dimensionalize time
t_stb = t_stb/Tpc;
t_cfd = t_cfd/Tpc;


%% PLOT VELOCITY DATA
% Color data
cg = [177 129 11]/255;
cmo = [255,155,26]/255;
cgr = [0,0,0]/255;
cssb = [124 166 192]/255;
cetb = [91 104 112]/255;
cfrt = [46 175 155]/255;
color_opts = [cg;cssb;cgr;cfrt;cmo;cetb];

% Set the colors for each modality
cmri = cgr;
cstb = cg;
ccfd = cssb;
%cstb_light = [238,232,170]/255;
cstb_light = [255,215,32]/255;
ccfd_light = [135,206,235]/255;
cmri_light = [50 50 50]/255;
% Set subplot positions
c2 = 0.50;%4846;
r2 = 0.2624;
% Set maximum flow rates for plotting
if strcmp(cadat.FILE.anynum,'111')
    max_plot_flow_1 = 9;% 0.5;
    inc_plot_flow_1 = 3;
    max_plot_flow_2 = 3;%0.18;
    inc_plot_flow_2 = 1;%0.06;
else
    max_plot_flow_1 = 9;
    inc_plot_flow_1 = 3;
    max_plot_flow_2 = 4;
    inc_plot_flow_2 = 1;
end
%%% PLOT FLOW RATE COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_plot_nums = [2,1,3,4];
close all

% Plot all inlets and outlets
for q = 1:1:size(q_mri,2)
    if q == 1
        figure(301); hold off; plot(t_mri,q_mri(:,q),'Color',cmri,'LineWidth',4);
        hold on; plot(t_stb,q_stb(:,q),'Color',cstb,'LineWidth',4);
        plot(t_cfd,q_cfd(:,q),'Color',ccfd,'LineWidth',4);
        set(gca,'FontSize',28); set(gcf,'Color',[1 1 1]);
        axis([0 1 0 0.4])
        set(gca,'XTick',0:0.2:1)
        set(gca,'YTick',0:inc_plot_flow_1:max_plot_flow_1)
        xlabel('t/T');
        ylabel('flow rate (mL/s)')
        pbaspect([1 0.6 1])
        h = legend('4D Flow','STB','CFD','Location','northoutside','Orientation','horizontal');
        set(h,'Box','off')
        print(301,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_comp_inlet.pdf'])
    else
        figure(302); hold off; h1 = subplot(2,2,sub_plot_nums(q-1)); plot(t_mri,q_mri(:,q),'Color',cmri,'LineWidth',2);
        hold on; plot(t_stb,q_stb(:,q),'Color',cstb,'LineWidth',2);
        plot(t_cfd,q_cfd(:,q),'Color',ccfd,'LineWidth',2);
        axis([0 1 0 max_plot_flow_2])
        set(gca,'FontSize',15); set(gcf,'Color',[1,1,1])
        pbaspect([1 0.6 1])
        % Set the y axis labels
        if sub_plot_nums(q-1) == 2 || sub_plot_nums(q-1) == 4
            % Turn off the y-axis
            set(gca,'YTickLabels',[]);
        else
            set(gca,'YTick',0:inc_plot_flow_2:max_plot_flow_2)
            %set(gca,'YTickLabels',0:0.03:0.18);
            ylabel('flow rate (mL/s)');
        end
        % Set the x axis labels
        if sub_plot_nums(q-1) == 2 || sub_plot_nums(q-1) == 1
            % Turn off the y-axis
            set(gca,'XTickLabels',[]);
        else
            set(gca,'XTick',0:0.2:1)
            %set(gca,'XTickLabels',0:0.03:0.18);
            xlabel('t/T');
        end
        % Set the position of the subplot
        if sub_plot_nums(q-1) == 2 || sub_plot_nums(q-1) == 4
            p1 = get(h1,'Position');
            p1(1) = c2;
            set(h1,'Position',p1)
        end
        if sub_plot_nums(q-1) == 3 || sub_plot_nums(q-1) == 4
            p1 = get(h1,'Position');
            p1(2) = r2;
            set(h1,'Position',p1)
        end
    end
end
% Print the flow rate plots
print(302,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_comp_outlet.pdf'])


% Compute flow rate in vs flow rate out
q_mri_in = q_mri(:,1);
q_mri_out = sum(q_mri(:,2:end),2);
q_mri_error = abs(q_mri_in - q_mri_out);
q_stb_in = q_stb(:,1);
q_stb_out = sum(q_stb(:,2:end),2);
q_stb_error = abs(q_stb_in - q_stb_out);
q_cfd_in = q_cfd(:,1);
q_cfd_out = sum(q_cfd(:,2:end),2);
q_cfd_error = abs(q_cfd_in - q_cfd_out);
        
% Plot flow rate using full resolution for CFD and STB
% Plot all inlets and outlets
line_size = 2;
font_size = 20;
c_base = 0.13; dc = 0.3;
r_base = 0.7093; dr = 0.28;
if strcmp(cadat.FILE.anynum,'111')
    y_step = [3,3,1,1,1];
    max_y_val = [9,9,3,3,3];
else
    y_step = [2,2,1,1];
    max_y_val = [8,8,4,4];
end

for q = 1:1:size(q_mri,2)
    if q == 1
        figure(401); hold off; plot(t_mri,q_mri(:,q),'Color',cmri,'LineWidth',4);
        hold on; plot(t_stb,q_full_stb(:,q),'Color',cstb_light,'LineWidth',4);
        plot(t_cfd,q_full_cfd(:,q),'Color',ccfd_light,'LineWidth',4);
        plot(t_stb,q_stb(:,q),'--','Color',cstb,'LineWidth',4);
        plot(t_cfd,q_cfd(:,q),'--','Color',ccfd,'LineWidth',4);
        set(gca,'FontSize',28); set(gcf,'Color',[1 1 1]);
        axis([0 1 0 max_plot_flow_1])
        set(gca,'XTick',0:0.2:1)
        set(gca,'YTick',0:inc_plot_flow_1:max_plot_flow_1)
        xlabel('t/T');
        ylabel('flow rate (mL/s)')
        pbaspect([1 0.6 1])
        h = legend('4D Flow','STB-Full','CFD-Full','Location','northoutside','Orientation','horizontal');
        set(h,'Box','off')
        print(401,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_FULLcomp_inlet.pdf'])
        savefig(401,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_FULLcomp_inlet.fig'])
        h = legend('','','','STB-Voxel','CFD-Voxel','Location','northoutside','Orientation','horizontal');
        set(h,'Box','off')
        print(401,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_FULLcomp_inlet_legend2.pdf'])
        plot([0.5265,0.5265]/Tpc,[0 max_plot_flow_1],'--','Color','r','LineWidth',2)
        plot([0.324,0.324]/Tpc,[0 max_plot_flow_1],'--','Color','r','LineWidth',2)
        plot([0.081,0.081]/Tpc,[0 max_plot_flow_1],'--','Color','r','LineWidth',2)
        h = legend('','','','','','','Location','northoutside','Orientation','horizontal');
        set(h,'Box','off');
        pbaspect([1 0.25 1])
        set(gca,'FontSize',20)
        print(401,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_time_points.pdf'])
    else
        figure(402); hold off; h1 = subplot(2,2,sub_plot_nums(q-1)); plot(t_mri,q_mri(:,q),'Color',cmri,'LineWidth',2);
        hold on; plot(t_stb,q_stb(:,q),'Color',cstb_light,'LineWidth',2);
        plot(t_cfd,q_cfd(:,q),'Color',ccfd_light,'LineWidth',2);
        plot(t_stb,q_full_stb(:,q),'--','Color',cstb,'LineWidth',2);
        plot(t_cfd,q_full_cfd(:,q),'--','Color',ccfd,'LineWidth',2);
        axis([0 1 0 max_plot_flow_2])
        set(gca,'FontSize',15); set(gcf,'Color',[1,1,1])
        pbaspect([1 0.6 1])
        % Set the y axis labels
        if sub_plot_nums(q-1) == 2 || sub_plot_nums(q-1) == 4
            % Turn off the y-axis
            set(gca,'YTickLabels',[]);
        else
            set(gca,'YTick',0:inc_plot_flow_2:max_plot_flow_2)
            %set(gca,'YTickLabels',0:0.03:0.18);
            ylabel('flow rate (mL/s)');
        end
        % Set the x axis labels
        if sub_plot_nums(q-1) == 2 || sub_plot_nums(q-1) == 1
            % Turn off the y-axis
            set(gca,'XTickLabels',[]);
        else
            set(gca,'XTick',0:0.2:1)
            %set(gca,'XTickLabels',0:0.03:0.18);
            xlabel('t/T');
        end
        % Set the position of the subplot
        if sub_plot_nums(q-1) == 2 || sub_plot_nums(q-1) == 4
            p1 = get(h1,'Position');
            p1(1) = c2;
            set(h1,'Position',p1)
        end
        if sub_plot_nums(q-1) == 3 || sub_plot_nums(q-1) == 4
            p1 = get(h1,'Position');
            p1(2) = r2;
            set(h1,'Position',p1)
        end
    end
    
    figure(403); hold off; h1 = subplot(3,2,q); 
    plot(t_mri,q_mri(:,q),'Color',cmri,'LineWidth',line_size);
    hold on; plot(t_stb,q_full_stb(:,q),'Color',cstb_light,'LineWidth',line_size);
    plot(t_cfd,q_full_cfd(:,q),'Color',ccfd_light,'LineWidth',line_size);
    plot(t_stb,q_stb(:,q),'--','Color',cstb,'LineWidth',line_size);
    plot(t_cfd,q_cfd(:,q),'--','Color',ccfd,'LineWidth',line_size);
    set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]);
    axis([0 1 0 max_y_val(q)])
    set(gca,'XTick',0:0.25:1)
    set(gca,'YTick',0:y_step(q):max_y_val(q))
    pbaspect([1 0.6 1]);
    p1 = get(h1,'Position');
    % Set the x-axis labels
    if q < 3
        set(gca,'XTickLabels',[]);
    end
    % Set the y-axis labels
    if (q == 2) || (q == 4) || (q == 6)
        set(gca,'YTickLabels',[]);
    end
    % Set the col position of subplots
    if (q == 2) || (q == 4) || (q == 6)
        p1(1) = c_base + dc;
    end
    % Set the y-axis labels
    if (q == 3) || (q == 4)
        p1(2) = r_base - dr;
    elseif (q == 5) || (q == 6)
        p1(2) = r_base - 2*dr;
    end
    set(h1,'Position',p1);
    
end
% Print the flow rate plots
print(402,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_FULLcomp_outlet.pdf'])
savefig(402,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_FULLcomp_outlet.fig'])

% Compute flow rate in vs flow rate out
q_full_stb_in = q_full_stb(:,1);
q_full_stb_out = sum(q_full_stb(:,2:end),2);
q_full_stb_error = abs(q_full_stb_in - q_full_stb_out);
q_full_cfd_in = q_full_cfd(:,1);
q_full_cfd_out = sum(q_full_cfd(:,2:end),2);
q_full_cfd_error = abs(q_full_cfd_in - q_full_cfd_out);

avg_cfd = mean(q_full_cfd,1);
avg_stbva = mean(q_stb,1);
avg_stb = mean(q_full_stb,1);
avg_mri = mean(q_mri,1);
avg_cfdva = mean(q_cfd,1);
avg_all = transpose([avg_mri;avg_stb;avg_cfd;avg_stbva;avg_cfdva]);

max_cfd = max(q_full_cfd,[],1);
max_stbva = max(q_stb,[],1);
max_stb = max(q_full_stb,[],1);
max_mri = max(q_mri,[],1);
max_cfdva = max(q_cfd,[],1);
max_all = transpose([max_mri;max_stb;max_cfd;max_stbva;max_cfdva]);

min_cfd = min(q_full_cfd,[],1);
min_stbva = min(q_stb,[],1);
min_stb = min(q_full_stb,[],1);
min_mri = min(q_mri,[],1);
min_cfdva = min(q_cfd,[],1);
min_all = transpose([min_mri;min_stb;min_cfd;min_stbva;min_cfdva]);

% Plot flow rate in vs flow rate out
line_size = 3;
figure(601); hold off; plot(t_mri,q_mri_in,'--','Color',cmri,'LineWidth',line_size)
hold on; plot(t_mri,q_mri_out,'-','Color',cmri_light,'LineWidth',line_size)
plot(t_stb,q_stb_in,'--','Color',cstb,'LineWidth',line_size)
hold on; plot(t_stb,q_stb_out,'-','Color',cstb_light,'LineWidth',line_size)
plot(t_cfd,q_cfd_in,'--','Color',ccfd,'LineWidth',line_size)
hold on; plot(t_cfd,q_cfd_out,'-','Color',ccfd_light,'LineWidth',line_size)
set(gca,'FontSize',28); set(gcf,'Color',[1 1 1]);
axis([0 1 0 max_plot_flow_1])
set(gca,'XTick',0:0.2:1)
set(gca,'YTick',0:0.1:max_plot_flow_1)
xlabel('t/T');
ylabel('flow rate (L/min)')
pbaspect([1 0.6 1])
h = legend('4D Flow-In','4D Flow-Out','STB-In','Location','northoutside','Orientation','horizontal');
set(h,'Box','off')
print(601,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_in_vs_out.pdf'])
h = legend('','','','STB-Out','CFD-In','CFD-Out','Location','northoutside','Orientation','horizontal');
set(h,'Box','off')
print(601,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_in_vs_out_legend2.pdf'])

% Plot flow rate error
figure(501); hold off; plot(t_mri,100*q_mri_error/max(q_mri(:)),'Color',cmri,'LineWidth',line_size);
hold on; plot(t_stb,100*q_full_stb_error/max(q_full_stb(:)),'Color',cstb_light,'LineWidth',line_size);
plot(t_cfd,100*q_full_cfd_error/max(q_full_cfd(:)),'Color',ccfd_light,'LineWidth',line_size);
plot(t_stb,100*q_stb_error/max(q_stb(:)),'--','Color',cstb,'LineWidth',line_size);
plot(t_cfd,100*q_cfd_error/max(q_cfd(:)),'--','Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',28); set(gcf,'Color',[1 1 1]);
axis([0 1 0 50])
set(gca,'XTick',0:0.2:1)
set(gca,'YTick',0:10:50)
xlabel('t/T');
ylabel('flow rate error (%)')
pbaspect([1 0.6 1])
h = legend('4D Flow','STB-Full','CFD-Full','Location','northoutside','Orientation','horizontal');
set(h,'Box','off')

print(501,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_error.pdf'])
savefig(501,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_error.fig'])


h = legend('','','','STB-Voxel','CFD-Voxel','Location','northoutside','Orientation','horizontal');
set(h,'Box','off')
print(501,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_flow_rate_error_legend2.pdf'])

max_flow_perc_error = [max(100*q_mri_error/max(q_mri(:))),max(100*q_full_stb_error/max(q_full_stb(:))),...
    max(100*q_full_cfd_error/max(q_full_cfd(:))),max(100*q_stb_error/max(q_stb(:))),max(100*q_cfd_error/max(q_cfd(:)))];
avg_flow_perc_error = [mean(100*q_mri_error/max(q_mri(:))),mean(100*q_full_stb_error/max(q_full_stb(:))),...
    mean(100*q_full_cfd_error/max(q_full_cfd(:))),mean(100*q_stb_error/max(q_stb(:))),mean(100*q_cfd_error/max(q_cfd(:)))];
max_flow_error = [max(q_mri_error),max(q_full_stb_error),max(q_full_cfd_error),max(q_stb_error),max(q_cfd_error)];
avg_flow_error = [mean(q_mri_error),mean(q_full_stb_error),mean(q_full_cfd_error),mean(q_stb_error),mean(q_cfd_error)];
flow_error_comp = [max_flow_error;max_flow_perc_error;avg_flow_error;avg_flow_perc_error];

%%% PLOT BLAND-ALTMAN COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plotting is done at peak systole, peak diastole, and mid diastole in 
% order to compare all modalities
%mean_color = [135,206,250]/255;
mean_color = [65,105,225]/255;
mean_color_2 = [65,161,225]/255;
% Iterate through the three cycle points
for q = 1:1:3
    % Iterate through all file types
    for ftype = 1:1:5
        % Import the velocity data for each image type
        if ftype == 1
            if strcmp(cadat.FILE.anynum,'111')
                veldir = '/Users/Melissa/Desktop/ANY/111/MRI/vel_registered_masked/';
                basename = 'any111_MRI_regmask_';
                curr_t = [3,9,14];
            else
                veldir = '/Users/Melissa/Desktop/ANY/007/MRI/vel_gridded/';
                basename = 'any007_MRI_regmask_gridded_';
                curr_t = [9,13,4];
            end
        elseif ftype == 2
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/STB/vel_p2s1_phaseavg/'];
            basename = ['any',cadat.FILE.anynum,'_STB_v1_p2s1_phaseavg_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [33,131,213];
            else
                curr_t = [240,359,102];
            end
        elseif ftype == 3
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/CFD/vel_gridded/'];
            basename = ['any',cadat.FILE.anynum,'_CFD_regmask_gridded_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [55,217,352];
            else
                curr_t = [240,359,102];
            end
        elseif ftype == 4
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/STB/vel_p3_voxavg/'];
            basename = ['any',cadat.FILE.anynum,'_STB_v1_p3_voxavg_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [33,131,213];
            else
                curr_t = [240,359,102];
            end
            
        else
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/CFD/vel_p3_voxavg/'];
            basename = ['any',cadat.FILE.anynum,'_CFD_p3_voxavg_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [55,217,352];
            else
                curr_t = [240,359,102];
            end
        end
        
        % Import the current time data
        fp = load([veldir,basename,num2str(curr_t(q),'%05i'),'.mat']);
        % Get the velocity field
        x = fp.x; y = fp.y; z = fp.z;
        u = fp.u; v = fp.v; w = fp.w;
        
        % Convert to column vectors if they are not
        if size(x,2)>1, x = x'; end
        if size(y,2)>1, y = y'; end
        if size(z,2)>1, z = z'; end
        if size(u,2)>1, u = u'; end
        if size(v,2)>1, v = v'; end
        if size(w,2)>1, w = w'; end
        
        % Save the data into the cell matrices
        x_all{ftype,1} = x;
        y_all{ftype,1} = y;
        z_all{ftype,1} = z;
        u_all{ftype,1} = u;
        v_all{ftype,1} = v;
        w_all{ftype,1} = w;
    end
    
    % Compare the data with all other data - the order of the B-A analysis
    % will be (1) MRI-STB Full res, (2) MRI-CFD Full res, (3) MRI-STB Voxel 
    % avg,(4) MRI-CFD Voxel avg, (5) STB Full res-CFD Full res, 
    % (6) STB Full res-STB Voxel avg, (7) STB Full res-CFD Voxel avg, 
    % (8) CFD Full res-STB Voxel avg, (9) CFD Full res-CFD Voxel avg, 
    % (10) STB Voxel avg-CFD Voxel avg
    if strcmp(cadat.FILE.anynum,'111')
        % Location comparison limits
        minX = -10.15; maxX = 7.55;
        minY = -5.15; maxY = 13.75;
        minZ = -7; maxZ = 7;
    else
        minX = -14; maxX = 8;
        minY = -6; maxY = 24;
        minZ = -10; maxZ = 15;
    end
    iter = 1; % Set the iteration counter
    % Subplot location settings
    dc = 0.19; c2 = 0.13+dc; c3 = 0.13+2*dc;
    dr = 0.25; r2 = 0.7093-dr; r3 = 0.7093-2*dr;
    %umin = -0.25; umax = 0.25; % Max and min for u-vel x-axis settings
    %wmin = -0.15; wmax = 0.15; % Max and min for u-vel x-axis settings
    %ymin = -0.25; ymax = 0.25; % Max and min for y-axis settings
    ymin = -0.44; ymax = 0.44; % Max and min for y-axis settings
    %umin = -0.4; umax = 0.4; % Max and min for u-vel x-axis settings
    for ftype1 = 1:1:2
        % Set the data for the first file type
        x1 = x_all(ftype1,1); x1 = x1{1,1};
        y1 = y_all(ftype1,1); y1 = y1{1,1};
        z1 = z_all(ftype1,1); z1 = z1{1,1};
        u1 = u_all(ftype1,1); u1 = u1{1,1};
        v1 = v_all(ftype1,1); v1 = v1{1,1};
        w1 = w_all(ftype1,1); w1 = w1{1,1};
        
        % Set the aneurysm mask points for the first time point
        if ftype1 == 1
            % Set the aneurysm mask for the first dataset
            any_mask_1 = mri_any_mask;
            % Set the velocity mask information for the current dataset
            xmask_1 = mri_maskX;
            ymask_1 = mri_maskY;
            zmask_1 = mri_maskZ;
        else
            any_mask_1 = stb_any_mask;
            % Set the velocity mask information for the current dataset
            xmask_1 = stb_maskX;
            ymask_1 = stb_maskY;
            zmask_1 = stb_maskZ;
        end
        
        % Iterate through all other modalities to compare
        for ftype2 = ftype1+1:1:3
            % Set the minimum threshold
            if ftype1 == 2 || ftype2 == 2 || ftype1 == 3 || ftype2 == 3
                grid_half = 0.15;
                dist_thresh = sqrt(3*grid_half^2);
            else
                grid_half = 0.67;
                dist_thresh = sqrt(3*grid_half^2);
            end
            
            % Set the aneurysm mask points for the first time point
            if ftype2 == 2
                % Set the aneurysm mask for the first dataset
                any_mask_2 = stb_any_mask;
                % Set the velocity mask information for the current dataset
                xmask_2 = stb_maskX;
                ymask_2 = stb_maskY;
                zmask_2 = stb_maskZ;
            else
                any_mask_2 = cfd_any_mask;
                % Set the velocity mask information for the current dataset
                xmask_2 = cfd_maskX;
                ymask_2 = cfd_maskY;
                zmask_2 = cfd_maskZ;
            end        
            
            % Set the data for the first file type
            x2 = x_all(ftype2,1); x2 = x2{1,1};
            y2 = y_all(ftype2,1); y2 = y2{1,1};
            z2 = z_all(ftype2,1); z2 = z2{1,1};
            u2 = u_all(ftype2,1); u2 = u2{1,1};
            v2 = v_all(ftype2,1); v2 = v2{1,1};
            w2 = w_all(ftype2,1); w2 = w2{1,1};
            
            % Set the primary array as the shorter one, secondary as longer
            % Set the aneurysm mask to be used based on the primary dataset
            if length(x1)>length(x2)
                xp = x2; yp = y2; zp = z2; up = u2; vp = v2; wp = w2; 
                xs = x1; ys = y1; zs = z1; us = u1; vs = v1; ws = w1;
                any_mask = any_mask_2;
                xmask = xmask_2; ymask = ymask_2; zmask = zmask_2;
                mult_factor = -1;
            else
                xs = x2; ys = y2; zs = z2; us = u2; vs = v2; ws = w2; 
                xp = x1; yp = y1; zp = z1; up = u1; vp = v1; wp = w1;
                any_mask = any_mask_1;
                xmask = xmask_1; ymask = ymask_1; zmask = zmask_1;
                mult_factor = 1;
            end
            
            %%% Iterate through all data points to create the B-A arrays
            ct = 1;
            % Initialize outputs
            dist_disparity = [];
            ba_vec1 = [];
            ba_vec2 = [];
            for pn = 1:1:length(xp)
                % Get the current point from the x primary
                cx = xp(pn);
                cy = yp(pn);
                cz = zp(pn);
                
                % Compute the distances from the current point
                pt_dists = sqrt((xs-cx).^2 + (ys-cy).^2 + (zs-cz).^2);
                [min_dist,sn] = min(pt_dists); % Get the minimum distance
                
                % Determine if the current value is in range
                in_range_x = (cx >= minX) & (cx <= maxX);
                in_range_y = (cy >= minY) & (cy <= maxY);
                in_range_z = (cz >= minZ) & (cz <= maxZ);
                all_in_range = in_range_x & in_range_y & in_range_z;
                
                % Determine if the current point is in the aneurysm sac
                % Find the closest point in the mask to the current point
                [~,min_ptX] = min(abs(xmask-cx)); % Find the closest mask point
                [~,min_ptY] = min(abs(ymask-cy)); % Find the closest mask point
                [~,min_ptZ] = min(abs(zmask-cz)); % Find the closest mask point
                xm_pt = xmask(min_ptX); 
                ym_pt = ymask(min_ptY);
                zm_pt = zmask(min_ptZ);
                % Determine if the current mask point is in the aneurysm
                % mask array
                m_pt = repmat([xm_pt,ym_pt,zm_pt],[size(any_mask,1),1]);
                m_pt_dist = sum(abs(m_pt-any_mask),2);
                min_pt_dist = min(m_pt_dist);
                in_any_mask = min_pt_dist < 0.25; % Determine if the point is in the aneurysm mask
                
                % If the minimum distance is below the threshold and the 
                % point is within the defined ROI, a match is accepted, add
                % the point to the output arrays
                if (min_dist < dist_thresh) && all_in_range && in_any_mask
                    dist_disparity(ct,1) = min_dist;
                    ba_vec1(ct,:) = [cx,cy,cz,up(pn),vp(pn),wp(pn)];
                    ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn),vs(sn),ws(sn)];
                    ct = ct + 1;
                end
            end
            
            % Set index skip for plotting
            if iter < 3, iskip = 1;
            else iskip = 10; end
            
            % Compute the B-A variables and  plot them
            pbrat = 1; % Size ratio of plots
            msize = 2; % Set the marker size
            line_size = 2; % Set the line width size
            font_size = 12; % Set the font size
            %%% vel-magnitude B-A analysis
            vmag1 = sqrt(ba_vec1(:,4).^2 + ba_vec1(:,5).^2 + ba_vec1(:,6).^2);
            vmag2 = sqrt(ba_vec2(:,4).^2 + ba_vec2(:,5).^2 + ba_vec2(:,6).^2);
            ba_mean = mean([vmag1,vmag2],2);
            ba_diff = vmag1 - vmag2;
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            %figure(200+iter); subplot(1,4,1)
            %hold off; plot(ba_mean,ba_diff,'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            %plot([-1 1], [CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            %plot([-1 1], [CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',[135,206,250]/255,'LineWidth',line_size); %%%plot the linear fit
            %axis([0 0.4 -0.5 0.5]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);

            %%% u-velocity B-A analysis
            ba_mean = mean([ba_vec1(:,4),ba_vec2(:,4)],2);
            ba_diff = mult_factor*(ba_vec1(:,4) - ba_vec2(:,4));
            %ba_diff = abs(ba_diff);
            %ba_diff = ba_diff_act./ba_mean;
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            figure(600+q); h1 = subplot(3,3,3*(iter-1)+1);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-1 0 1], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-1 0 1], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-1 0 1], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
            axis([-0.44 0.44 ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'XTick',-0.4:0.2:0.4); grid('on');
            ylabel('difference (m/s)')
            % Set the axis information
            if iter < 3
                set(gca,'XTickLabel',{});
            else
                xlabel('mean (m/s)')
            end
            p = get(h1,'Position');
            % Set the subplot location (row)
            if 3*(iter-1)+1 > 6, p(2) = r3;
            elseif 3*(iter-1)+1 > 3, p(2) = r2;
            end
            set(h1,'Position',p);
            % Set title
            if iter == 1, title(h1,'u-velocity'); end
            % Save info in array
            ba_results{q,1}(:,3*(iter-1)+1) = [meanDiff;CR(1);CR(2)];
            
            %%% v-velocity B-A analysis
            ba_mean = mean([ba_vec1(:,5),ba_vec2(:,5)],2);
            ba_diff = mult_factor*(ba_vec1(:,5) - ba_vec2(:,5));
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            figure(600+q); h1 = subplot(3,3,3*(iter-1)+2);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-1 0 1], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-1 0 1], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-1 0 1], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
            axis([-0.25 0.45 ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'YTickLabel',{}); set(gca,'XTick',-0.2:0.2:0.4); grid('on');
            % Set the axis information
            if iter < 3
                set(gca,'XTickLabel',{});
            else
                xlabel('mean (m/s)')
            end
            p = get(h1,'Position');
            % Set the subplot location (row)
            if 3*(iter-1)+2 > 6, p(2) = r3;
            elseif 3*(iter-1)+2 > 3, p(2) = r2;
            end
            p(1) = c2; % Set column
            set(h1,'Position',p);
            if iter == 1, title(h1,'v-velocity'); end
            % Save info in array
            ba_results{q,1}(:,3*(iter-1)+2) = [meanDiff;CR(1);CR(2)];
            
            %%% w-velocity B-A analysis
            ba_mean = mean([ba_vec1(:,6),ba_vec2(:,6)],2);
            ba_diff = mult_factor*(ba_vec1(:,6) - ba_vec2(:,6));
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            figure(600+q); h1 = subplot(3,3,3*(iter-1)+3);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-1 0 1], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-1 0 1], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-1 0 1], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
            axis([-0.25 0.25 ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'YTickLabel',{}); set(gca,'XTick',-0.2:0.2:0.2); grid('on');
            % Set the axis information
            if iter < 3
                set(gca,'XTickLabel',{});
            else
                xlabel('mean (m/s)')
            end
            p = get(h1,'Position');
            % Set the subplot location (row)
            if 3*(iter-1)+3 > 6, p(2) = r3;
            elseif 3*(iter-1)+3 > 3, p(2) = r2;
            end
            p(1) = c3; % Set column
            set(h1,'Position',p);
            if iter == 1, title(h1,'w-velocity'); end
            % Save info in array
            ba_results{q,1}(:,3*(iter-1)+3) = [meanDiff;CR(1);CR(2)];
            
            % Save the disparity error in an array
            dist_disparity_error{iter,1} = dist_disparity;
            dist_disparity_avgerror(iter,1) = mean(dist_disparity);        
            
            % Increment the iteration counter
            iter = iter + 1;    
        end
    end
    % Get the cycle point print info
    if q == 1, cycle_pt = 'MidDiastole_';
    elseif q == 2, cycle_pt = 'PeakDiastole_';
    else, cycle_pt = 'PeakSystole_'; 
    end   
    % Print the Bland-Altman plots
    print(600+q,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-anymask-comp.pdf'])
end
% Print the legend
figure(600+q+1);
hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
plot([-1 1], [CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
plot([-1 0 1], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
axis([-0.25 0.25 ymin ymax]); set(gca,'FontSize',18); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
set(gca,'YTickLabel',{}); set(gca,'XTick',-0.2:0.2:0.2); grid('on');
h = legend('Data Points','95% CI','Linear Fit','Mean Diff','Location','northoutside','Orientation','horizontal');
set(h,'Box','off');
print(600+q+1,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_legend.pdf'])



%%% PLOT BLAND-ALTMAN ALL VELOCITY COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plotting is done at peak systole, peak diastole, and mid diastole in 
% order to compare all modalities. This is done using all velocity points
% as opposed to the previous analysis which is only in the aneurysm sac
%mean_color = [135,206,250]/255;
mean_color = [65,105,225]/255;
% Iterate through the three cycle points
for q = 3:1:3
    % Iterate through all file types
    for ftype = 1:1:5
        % Import the velocity data for each image type
        % Import the velocity data for each image type
        if ftype == 1
            if strcmp(cadat.FILE.anynum,'111')
                veldir = '/Users/Melissa/Desktop/ANY/111/MRI/vel_registered_masked/';
                basename = 'any111_MRI_regmask_';
                curr_t = [3,9,14];
            else
                veldir = '/Users/Melissa/Desktop/ANY/007/MRI/vel_gridded/';
                basename = 'any007_MRI_regmask_gridded_';
                curr_t = [9,13,4];
            end
        elseif ftype == 2
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/STB/vel_p2s1_phaseavg/'];
            basename = ['any',cadat.FILE.anynum,'_STB_v1_p2s1_phaseavg_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [33,131,213];
            else
                curr_t = [240,359,102];
            end
        elseif ftype == 3
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/CFD/vel_gridded/'];
            basename = ['any',cadat.FILE.anynum,'_CFD_regmask_gridded_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [55,217,352];
            else
                curr_t = [240,359,102];
            end
        elseif ftype == 4
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/STB/vel_p3_voxavg/'];
            basename = ['any',cadat.FILE.anynum,'_STB_v1_p3_voxavg_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [33,131,213];
            else
                curr_t = [240,359,102];
            end
            
        else
            veldir = ['/Users/Melissa/Desktop/ANY/',cadat.FILE.anynum,'/CFD/vel_p3_voxavg/'];
            basename = ['any',cadat.FILE.anynum,'_CFD_p3_voxavg_'];
            if strcmp(cadat.FILE.anynum,'111')
                curr_t = [55,217,352];
            else
                curr_t = [240,359,102];
            end
        end
        
        % Import the current time data
        fp = load([veldir,basename,num2str(curr_t(q),'%05i'),'.mat']);
        % Get the velocity field
        x = fp.x; y = fp.y; z = fp.z;
        u = fp.u; v = fp.v; w = fp.w;
        
        % Convert to column vectors if they are not
        if size(x,2)>1, x = x'; end
        if size(y,2)>1, y = y'; end
        if size(z,2)>1, z = z'; end
        if size(u,2)>1, u = u'; end
        if size(v,2)>1, v = v'; end
        if size(w,2)>1, w = w'; end
        
        % Save the data into the cell matrices
        x_all{ftype,1} = x;
        y_all{ftype,1} = y;
        z_all{ftype,1} = z;
        u_all{ftype,1} = u;
        v_all{ftype,1} = v;
        w_all{ftype,1} = w;
    end
    
    % Compare the data with all other data - the order of the B-A analysis
    % will be (1) MRI-STB Full res, (2) MRI-CFD Full res, (3) MRI-STB Voxel 
    % avg,(4) MRI-CFD Voxel avg, (5) STB Full res-CFD Full res, 
    % (6) STB Full res-STB Voxel avg, (7) STB Full res-CFD Voxel avg, 
    % (8) CFD Full res-STB Voxel avg, (9) CFD Full res-CFD Voxel avg, 
    % (10) STB Voxel avg-CFD Voxel avg
    if strcmp(cadat.FILE.anynum,'111')
        % Location comparison limits
        minX = -10.15; maxX = 7.55;
        minY = -5.15; maxY = 13.75;
        minZ = -7; maxZ = 7;
    else
        minX = -14; maxX = 8;
        minY = -6; maxY = 24;
        minZ = -10; maxZ = 15;
    end
    iter = 1; % Set the iteration counter
    % Subplot location settings
    dc = 0.19; c2 = 0.13+dc; c3 = 0.13+2*dc;
    dr = 0.25; r2 = 0.7093-dr; r3 = 0.7093-2*dr;
    %umin = -0.25; umax = 0.25; % Max and min for u-vel x-axis settings
    %wmin = -0.15; wmax = 0.15; % Max and min for u-vel x-axis settings
    %ymin = -0.25; ymax = 0.25; % Max and min for y-axis settings
    ymin = -0.44; ymax = 0.44; % Max and min for y-axis settings
    %umin = -0.4; umax = 0.4; % Max and min for u-vel x-axis settings
    for ftype1 = 1:1:2
        % Set the data for the first file type
        x1 = x_all(ftype1,1); x1 = x1{1,1};
        y1 = y_all(ftype1,1); y1 = y1{1,1};
        z1 = z_all(ftype1,1); z1 = z1{1,1};
        u1 = u_all(ftype1,1); u1 = u1{1,1};
        v1 = v_all(ftype1,1); v1 = v1{1,1};
        w1 = w_all(ftype1,1); w1 = w1{1,1};
        % Iterate through all other modalities to compare
        for ftype2 = ftype1+1:1:3
            % Set the minimum threshold
            if ftype1 == 2 || ftype2 == 2 || ftype1 == 3 || ftype2 == 3
                grid_half = 0.15;
                dist_thresh = sqrt(3*grid_half^2);
            else
                grid_half = 0.67;
                dist_thresh = sqrt(3*grid_half^2);
            end
            
            % Set the data for the first file type
            x2 = x_all(ftype2,1); x2 = x2{1,1};
            y2 = y_all(ftype2,1); y2 = y2{1,1};
            z2 = z_all(ftype2,1); z2 = z2{1,1};
            u2 = u_all(ftype2,1); u2 = u2{1,1};
            v2 = v_all(ftype2,1); v2 = v2{1,1};
            w2 = w_all(ftype2,1); w2 = w2{1,1};
            
            % Set the primary array as the shorter one, secondary as longer
            if length(x1)>length(x2)
                xp = x2; yp = y2; zp = z2; up = u2; vp = v2; wp = w2; 
                xs = x1; ys = y1; zs = z1; us = u1; vs = v1; ws = w1;
            else
                xs = x2; ys = y2; zs = z2; us = u2; vs = v2; ws = w2; 
                xp = x1; yp = y1; zp = z1; up = u1; vp = v1; wp = w1;
            end
            
            %%% Iterate through all data points to create the B-A arrays
            ct = 1;
            % Initialize outputs
            dist_disparity = [];
            ba_vec1 = [];
            ba_vec2 = [];
            for pn = 1:1:length(xp)
                % Get the current point from the x primary
                cx = xp(pn);
                cy = yp(pn);
                cz = zp(pn);
                
                % Compute the distances from the current point
                pt_dists = sqrt((xs-cx).^2 + (ys-cy).^2 + (zs-cz).^2);
                [min_dist,sn] = min(pt_dists); % Get the minimum distance
                
                % Determine if the current value is in range
                in_range_x = (cx >= minX) & (cx <= maxX);
                in_range_y = (cy >= minY) & (cy <= maxY);
                in_range_z = (cz >= minZ) & (cz <= maxZ);
                all_in_range = in_range_x & in_range_y & in_range_z;
                % If the minimum distance is below the threshold and the 
                % point is within the defined ROI, a match is accepted, add
                % the point to the output arrays
                if (min_dist < dist_thresh) && all_in_range
                    dist_disparity(ct,1) = min_dist;
                    ba_vec1(ct,:) = [cx,cy,cz,up(pn),vp(pn),wp(pn)];
                    ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn),vs(sn),ws(sn)];
                    ct = ct + 1;
                end
            end
            
            % Set index skip for plotting
            if iter < 3, iskip = 1;
            else iskip = 10; end
            
            % Compute the B-A variables and  plot them
            pbrat = 1; % Size ratio of plots
            msize = 2; % Set the marker size
            line_size = 2; % Set the line width size
            font_size = 12; % Set the font size
            %%% vel-magnitude B-A analysis
            vmag1 = sqrt(ba_vec1(:,4).^2 + ba_vec1(:,5).^2 + ba_vec1(:,6).^2);
            vmag2 = sqrt(ba_vec2(:,4).^2 + ba_vec2(:,5).^2 + ba_vec2(:,6).^2);
            ba_mean = mean([vmag1,vmag2],2);
            ba_diff = vmag1 - vmag2;
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            %figure(200+iter); subplot(1,4,1)
            %hold off; plot(ba_mean,ba_diff,'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            %plot([-1 1], [CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            %plot([-1 1], [CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',[135,206,250]/255,'LineWidth',line_size); %%%plot the linear fit
            %axis([0 0.4 -0.5 0.5]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);

            %%% u-velocity B-A analysis
            ba_mean = mean([ba_vec1(:,4),ba_vec2(:,4)],2);
            ba_diff = ba_vec1(:,4) - ba_vec2(:,4);
            %ba_diff = abs(ba_diff);
            %ba_diff = ba_diff_act./ba_mean;
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            figure(200+q); h1 = subplot(3,3,3*(iter-1)+1);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-1 0 1], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-1 0 1], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-1 0 1], [meanDiff meanDiff meanDiff],'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            axis([-0.44 0.44 ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'XTick',-0.4:0.2:0.4); grid('on');
            ylabel('difference (m/s)')
            % Set the axis information
            if iter < 3
                set(gca,'XTickLabel',{});
            else
                xlabel('mean (m/s)')
            end
            p = get(h1,'Position');
            % Set the subplot location (row)
            if 3*(iter-1)+1 > 6, p(2) = r3;
            elseif 3*(iter-1)+1 > 3, p(2) = r2;
            end
            set(h1,'Position',p);
            % Set title
            if iter == 1, title('u-velocity'); end
            
            %%% v-velocity B-A analysis
            ba_mean = mean([ba_vec1(:,5),ba_vec2(:,5)],2);
            ba_diff = ba_vec1(:,5) - ba_vec2(:,5);
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            figure(200+q); h1 = subplot(3,3,3*(iter-1)+2);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-1 0 1], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-1 0 1], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-1 0 1], [meanDiff meanDiff meanDiff],'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            axis([-0.25 0.45 ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'YTickLabel',{}); set(gca,'XTick',-0.2:0.2:0.4); grid('on');
            % Set the axis information
            if iter < 3
                set(gca,'XTickLabel',{});
            else
                xlabel('mean (m/s)')
            end
            p = get(h1,'Position');
            % Set the subplot location (row)
            if 3*(iter-1)+2 > 6, p(2) = r3;
            elseif 3*(iter-1)+2 > 3, p(2) = r2;
            end
            p(1) = c2; % Set column
            set(h1,'Position',p);
            if iter == 1, title('v-velocity'); end
            
            %%% w-velocity B-A analysis
            ba_mean = mean([ba_vec1(:,6),ba_vec2(:,6)],2);
            ba_diff = ba_vec1(:,6) - ba_vec2(:,6);
            meanDiff = mean(ba_diff);
            sdDiff = std(ba_diff);
            CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
            % Compute linear fit of means and differences
            linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
            % Plot the u-velocity B-A analysis
            figure(200+q); h1 = subplot(3,3,3*(iter-1)+3);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-1 0 1], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-1 0 1], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-1 0 1], [meanDiff meanDiff meanDiff],'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            axis([-0.25 0.25 ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'YTickLabel',{}); set(gca,'XTick',-0.2:0.2:0.2); grid('on');
            % Set the axis information
            if iter < 3
                set(gca,'XTickLabel',{});
            else
                xlabel('mean (m/s)')
            end
            p = get(h1,'Position');
            % Set the subplot location (row)
            if 3*(iter-1)+3 > 6, p(2) = r3;
            elseif 3*(iter-1)+3 > 3, p(2) = r2;
            end
            p(1) = c3; % Set column
            set(h1,'Position',p);
            if iter == 1, title('w-velocity'); end
            
            % Save the disparity error in an array
            dist_disparity_error{iter,1} = dist_disparity;
            dist_disparity_avgerror(iter,1) = mean(dist_disparity);        
            
            % Increment the iteration counter
            iter = iter + 1;    
        end
    end
    % Get the cycle point print info
    if q == 1, cycle_pt = 'MidDiastole_';
    elseif q == 2, cycle_pt = 'PeakDiastole_';
    else, cycle_pt = 'PeakSystole_'; 
    end   
    % Print the Bland-Altman plots
    print(200+q,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-comp.pdf'])
end
% Print the legend
%figure(200+q+1);
%hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
%plot([-1 1], [CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
%plot([-1 1], [-1 1].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
%axis([-0.25 0.25 ymin ymax]); set(gca,'FontSize',18); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
%set(gca,'YTickLabel',{}); set(gca,'XTick',-0.2:0.2:0.2); grid('on');
%h = legend('Data Points','95% CI','Mean Diff','Location','northoutside','Orientation','horizontal');
%set(h,'Box','off');
%print(200+q+1,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_legend.pdf'])



%%% PLOT VELOCITY DISTRIBUTION COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The plotting here is done for peak systole and for all velocities across
% the entire cycle
% Iterate through all mask points and extract the current 

% MRI mask velocity points
ct = 1;
for t = mri_file_s:1:mri_file_e
    % Load MRI velocity field
    if strcmp(cadat.FILE.anynum,'007')
        fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'MRI',fslash,'vel_gridded',fslash,'any',cadat.FILE.anynum,'_MRI_regmask_gridded_',num2str(t,'%05i'),'.mat']);
    else
        fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'MRI',fslash,'vel_registered_masked',fslash,'any',cadat.FILE.anynum,'_MRI_regmask_',num2str(t,'%05i'),'.mat']);
    end
    vel_x = fp.x; vel_y = fp.y; vel_z = fp.z;
    vel_u = fp.u; vel_v = fp.v; vel_w = fp.w;
    for q = 1:1:size(mri_any_mask,1)
        % Get the current x, y, z points in the mask
        cx = mri_any_mask(q,1);
        cy = mri_any_mask(q,2);
        cz = mri_any_mask(q,3);
        % Find the point in the current field closest to this
        pt_dist = sqrt((vel_x-cx).^2 + (vel_y-cy).^2 + (vel_z-cz).^2);
        [~,close_ind] = min(pt_dist);
        mri_u_dis(ct,1) = vel_u(close_ind);
        mri_v_dis(ct,1) = vel_v(close_ind);
        mri_w_dis(ct,1) = vel_w(close_ind);
        ct = ct + 1;
    end
end
mri_vel_dis = sqrt(mri_u_dis.^2 + mri_v_dis.^2 + mri_w_dis.^2);

% STB mask velocity points
ct = 1;
for t = stb_file_s:stb_file_skip:stb_file_e
    % Load STB velocity field
    fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'STB',fslash,'vel_p3_voxavg',fslash,'any',cadat.FILE.anynum,'_STB_v1_p3_voxavg_',num2str(t,'%05i'),'.mat']);
    vel_x = fp.x; vel_y = fp.y; vel_z = fp.z;
    vel_u = fp.u; vel_v = fp.v; vel_w = fp.w;
    for q = 1:1:size(stb_any_mask,1)
        % Get the current x, y, z points in the mask
        cx = stb_any_mask(q,1);
        cy = stb_any_mask(q,2);
        cz = stb_any_mask(q,3);
        % Find the point in the current field closest to this
        pt_dist = sqrt((vel_x-cx).^2 + (vel_y-cy).^2 + (vel_z-cz).^2);
        [~,close_ind] = min(pt_dist);
        stb_u_dis(ct,1) = vel_u(close_ind);
        stb_v_dis(ct,1) = vel_v(close_ind);
        stb_w_dis(ct,1) = vel_w(close_ind);
        ct = ct + 1;
    end
end
stb_vel_dis = sqrt(stb_u_dis.^2 + stb_v_dis.^2 + stb_w_dis.^2);

% CFD mask velocity points
ct = 1;
for t = cfd_file_s:1:cfd_file_e
    % Load STB velocity field
    fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'CFD',fslash,'vel_p3_voxavg',fslash,'any',cadat.FILE.anynum,'_CFD_p3_voxavg_',num2str(t,'%05i'),'.mat']);
    vel_x = fp.x; vel_y = fp.y; vel_z = fp.z;
    vel_u = fp.u; vel_v = fp.v; vel_w = fp.w;
    for q = 1:1:size(cfd_any_mask,1)
        % Get the current x, y, z points in the mask
        cx = cfd_any_mask(q,1);
        cy = cfd_any_mask(q,2);
        cz = cfd_any_mask(q,3);
        % Find the point in the current field closest to this
        pt_dist = sqrt((vel_x-cx).^2 + (vel_y-cy).^2 + (vel_z-cz).^2);
        [~,close_ind] = min(pt_dist);
        cfd_u_dis(ct,1) = vel_u(close_ind);
        cfd_v_dis(ct,1) = vel_v(close_ind);
        cfd_w_dis(ct,1) = vel_w(close_ind);
        ct = ct + 1;
    end
end
cfd_vel_dis = sqrt(cfd_u_dis.^2 + cfd_v_dis.^2 + cfd_w_dis.^2);


if CREATE_FRES_DIS
    % STB mask velocity points (STB full resolution)
    % Using the velocity mask, identify all points who's closest point on the
    % voxel averaged field is in the aneurysm mask
    % Start by identifying the voxel averaged mask points that are in the
    % aneurysm
    [Xva,Yva,Zva] = meshgrid(stb_maskX,stb_maskY,stb_maskZ);
    stb_vamask = reshape(stb_velmask,[size(stb_velmask,1)*size(stb_velmask,2)*size(stb_velmask,3),1]);
    stb_vamaskX = reshape(Xva,[size(stb_velmask,1)*size(stb_velmask,2)*size(stb_velmask,3),1]);
    stb_vamaskY = reshape(Yva,[size(stb_velmask,1)*size(stb_velmask,2)*size(stb_velmask,3),1]);
    stb_vamaskZ = reshape(Zva,[size(stb_velmask,1)*size(stb_velmask,2)*size(stb_velmask,3),1]);
    stb_vamask_inany = zeros(size(stb_vamask));
    for q = 1:1:size(stb_any_mask,1)
        % Determine which points in the va mask the aneurysm points correspond
        % to
        cx = stb_any_mask(q,1);
        cy = stb_any_mask(q,2);
        cz = stb_any_mask(q,3);
        % Find the closest point on the va mask
        pt_dist = sqrt((stb_vamaskX-cx).^2 + (stb_vamaskY-cy).^2 + (stb_vamaskZ-cz).^2);
        [~,min_ind] = min(pt_dist);
        % Set the point as in the mask
        stb_vamask_inany(min_ind) = 1;
    end
    % Iterate through all time
    stb_fres_u_dis = [];
    stb_fres_v_dis = [];
    stb_fres_w_dis = [];
    for t = stb_file_s:stb_file_skip:stb_file_e
        % Load the full resolution STB file
        fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'STB',fslash,'vel_p2s1_phaseavg',fslash,'any',cadat.FILE.anynum,'_STB_v1_p2s1_phaseavg_',num2str(t,'%05i'),'.mat']);
        vel_x = fp.x; vel_y = fp.y; vel_z = fp.z;
        vel_u = fp.u; vel_v = fp.v; vel_w = fp.w;
        % Iterate through each point in the full resolution wss vector. Add full
        % resolution point if its closest voxel averaged point is in the aneurysm
        % mask
        fres_inany = zeros(size(vel_x));
        for q = 1:1:size(vel_x)
            % Get current x,y,z point
            cx = vel_x(q);
            cy = vel_y(q);
            cz = vel_z(q);
            % Find the closest point on the va mask to the current point
            % Find the closest point on the va mask
            pt_dist = sqrt((stb_vamaskX-cx).^2 + (stb_vamaskY-cy).^2 + (stb_vamaskZ-cz).^2);
            [~,min_ind] = min(pt_dist);
            % Determine if the closest point is in the aneurysm mask
            in_mask = stb_vamask_inany(min_ind);
            % If it is in the mask, add it to the full resolution in mask
            if in_mask, fres_inany(q) = 1; end
        end
        % Keep only the wall shear stress values in the anuerysm
        stb_fres_u_dis = [stb_fres_u_dis;vel_u(fres_inany==1)];
        stb_fres_v_dis = [stb_fres_v_dis;vel_v(fres_inany==1)];
        stb_fres_w_dis = [stb_fres_w_dis;vel_w(fres_inany==1)];
    end
    stb_fres_vel_dis = sqrt(stb_fres_u_dis.^2 + stb_fres_v_dis.^2 + stb_fres_w_dis.^2);

    % CFD mask velocity points (CFD full resolution)
    % Using the velocity mask, identify all points who's closest point on the
    % voxel averaged field is in the aneurysm mask
    % Start by identifying the voxel averaged mask points that are in the
    % aneurysm
    [Xva,Yva,Zva] = meshgrid(cfd_maskX,cfd_maskY,cfd_maskZ);
    cfd_vamask = reshape(cfd_velmask,[size(cfd_velmask,1)*size(cfd_velmask,2)*size(cfd_velmask,3),1]);
    cfd_vamaskX = reshape(Xva,[size(cfd_velmask,1)*size(cfd_velmask,2)*size(cfd_velmask,3),1]);
    cfd_vamaskY = reshape(Yva,[size(cfd_velmask,1)*size(cfd_velmask,2)*size(cfd_velmask,3),1]);
    cfd_vamaskZ = reshape(Zva,[size(cfd_velmask,1)*size(cfd_velmask,2)*size(cfd_velmask,3),1]);
    cfd_vamask_inany = zeros(size(cfd_vamask));
    for q = 1:1:size(cfd_any_mask,1)
        % Determine which points in the va mask the aneurysm points correspond
        % to
        cx = cfd_any_mask(q,1);
        cy = cfd_any_mask(q,2);
        cz = cfd_any_mask(q,3);
        % Find the closest point on the va mask
        pt_dist = sqrt((cfd_vamaskX-cx).^2 + (cfd_vamaskY-cy).^2 + (cfd_vamaskZ-cz).^2);
        [~,min_ind] = min(pt_dist);
        % Set the point as in the mask
        cfd_vamask_inany(min_ind) = 1;
    end
    % Iterate through all time
    cfd_fres_u_dis = [];
    cfd_fres_v_dis = [];
    cfd_fres_w_dis = [];
    for t = cfd_file_s:cfd_file_skip:cfd_file_e
        % Load the full resolution STB file
        fp = load([cadat.DIR.dirbase,fslash,cadat.FILE.anynum,fslash,'CFD',fslash,'vel_gridded',fslash,'any',cadat.FILE.anynum,'_CFD_regmask_gridded_',num2str(t,'%05i'),'.mat']);
        vel_x = fp.x; vel_y = fp.y; vel_z = fp.z;
        vel_u = fp.u; vel_v = fp.v; vel_w = fp.w;
        % Iterate through each point in the full resolution wss vector. Add full
        % resolution point if its closest voxel averaged point is in the aneurysm
        % mask
        fres_inany = zeros(size(vel_x));
        for q = 1:1:size(vel_x)
            % Get current x,y,z point
            cx = vel_x(q);
            cy = vel_y(q);
            cz = vel_z(q);
            % Find the closest point on the va mask to the current point
            % Find the closest point on the va mask
            pt_dist = sqrt((cfd_vamaskX-cx).^2 + (cfd_vamaskY-cy).^2 + (cfd_vamaskZ-cz).^2);
            [~,min_ind] = min(pt_dist);
            % Determine if the closest point is in the aneurysm mask
            in_mask = cfd_vamask_inany(min_ind);
            % If it is in the mask, add it to the full resolution in mask
            if in_mask, fres_inany(q) = 1; end
        end
        % Keep only the wall shear stress values in the anuerysm
        cfd_fres_u_dis = [cfd_fres_u_dis;vel_u(fres_inany==1)];
        cfd_fres_v_dis = [cfd_fres_v_dis;vel_v(fres_inany==1)];
        cfd_fres_w_dis = [cfd_fres_w_dis;vel_w(fres_inany==1)];
    end
    cfd_fres_vel_dis = sqrt(cfd_fres_u_dis.^2 + cfd_fres_v_dis.^2 + cfd_fres_w_dis.^2);

    % Save the full resolution velocity distributions
    save([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'STB_fres_vel_dist','.mat'],'stb_fres_u_dis','stb_fres_v_dis','stb_fres_w_dis','stb_fres_vel_dis')
    save([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'CFD_fres_vel_dist','.mat'],'cfd_fres_u_dis','cfd_fres_v_dis','cfd_fres_w_dis','cfd_fres_vel_dis')
else
    % If they are not to be computed, load them
    fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'STB_fres_vel_dist','.mat']);
    stb_fres_u_dis = fp.stb_fres_u_dis;
    stb_fres_v_dis = fp.stb_fres_v_dis;
    stb_fres_w_dis = fp.stb_fres_w_dis;
    stb_fres_vel_dis = fp.stb_fres_vel_dis;
    fp = load([cadat.DIR.saveppfiles,'any',cadat.FILE.anynum,'CFD_fres_vel_dist','.mat']);
    cfd_fres_u_dis = fp.cfd_fres_u_dis;
    cfd_fres_v_dis = fp.cfd_fres_v_dis;
    cfd_fres_w_dis = fp.cfd_fres_w_dis;
    cfd_fres_vel_dis = fp.cfd_fres_vel_dis;
end

% Plot the velocity distribution
num_bins = 75;
fsize = 28;
line_size = 3;


% For the velocity magnitude
min_bin = min([min(cfd_vel_dis),min(stb_vel_dis),min(mri_vel_dis),min(cfd_vel_dis),min(stb_vel_dis)]);
max_bin = max([max(cfd_vel_dis),max(stb_vel_dis),max(mri_vel_dis),max(cfd_fres_vel_dis),max(stb_fres_vel_dis)]);
bin_edges = linspace(min_bin,max_bin,num_bins+1);
figure(1);
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
mri_pdf = histogram(mri_vel_dis,bin_edges,'Normalization','pdf');
mri_pdf = mri_pdf.Values;
mri_cdf = histogram(mri_vel_dis,bin_edges,'Normalization','cdf');
mri_cdf = mri_cdf.Values;
stb_pdf = histogram(stb_vel_dis,bin_edges,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_vel_dis,bin_edges,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_fres_vel_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_fres_pdf = histogram(cfd_fres_vel_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
cfd_pdf = histogram(cfd_vel_dis,bin_edges,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_vel_dis,bin_edges,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
max_val = max([mri_pdf,stb_pdf,cfd_pdf,stb_fres_pdf,cfd_fres_pdf]);
figure(305); hold off; plot(plt_bins,mri_pdf/max_val,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',line_size);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',line_size);
plot(plt_bins,stb_pdf/max_val,'--','Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
xlabel('velocity magnitude (m/s)')
ylabel('PDF')
axis([0 0.6 0 1])
figure(306); hold off; plot(plt_bins,mri_cdf,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_cdf,'Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_cdf,'Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
xlabel('velocity magnitude (m/s)')
ylabel('CDF')

% For the u velocity
min_bin = min([min(cfd_u_dis),min(stb_u_dis),min(mri_u_dis)]);
max_bin = max([max(cfd_u_dis),max(stb_u_dis),max(mri_u_dis)]);
bin_edges = linspace(min_bin,max_bin,num_bins+1);
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
figure(1);
mri_pdf = histogram(mri_u_dis,bin_edges,'Normalization','pdf');
mri_pdf = mri_pdf.Values;
mri_cdf = histogram(mri_u_dis,bin_edges,'Normalization','cdf');
mri_cdf = mri_cdf.Values;
stb_pdf = histogram(stb_u_dis,bin_edges,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_u_dis,bin_edges,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_fres_u_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_fres_pdf = histogram(cfd_fres_u_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
cfd_pdf = histogram(cfd_u_dis,bin_edges,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_u_dis,bin_edges,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
max_val = max([mri_pdf,stb_pdf,cfd_pdf,stb_fres_pdf,cfd_fres_pdf]);
figure(307); hold off; plot(plt_bins,mri_pdf/max_val,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',line_size);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',line_size);
plot(plt_bins,stb_pdf/max_val,'--','Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
xlabel('u velocity (m/s)')
axis([-0.3 0.3 0 1])
ylabel('PDF')
figure(308); hold off; plot(plt_bins,mri_cdf,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_cdf,'Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_cdf,'Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
xlabel('u velocity (m/s)')
ylabel('CDF')

% For v velocity magnitude
min_bin = min([min(cfd_v_dis),min(stb_v_dis),min(mri_v_dis)]);
max_bin = max([max(cfd_v_dis),max(stb_v_dis),max(mri_v_dis)]);
bin_edges = linspace(min_bin,max_bin,num_bins+1);
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
figure(1);
mri_pdf = histogram(mri_v_dis,bin_edges,'Normalization','pdf');
mri_pdf = mri_pdf.Values;
mri_cdf = histogram(mri_v_dis,bin_edges,'Normalization','cdf');
mri_cdf = mri_cdf.Values;
stb_pdf = histogram(stb_v_dis,bin_edges,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_v_dis,bin_edges,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_fres_v_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_fres_pdf = histogram(cfd_fres_v_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
cfd_pdf = histogram(cfd_v_dis,bin_edges,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_v_dis,bin_edges,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
max_val = max([mri_pdf,stb_pdf,cfd_pdf,stb_fres_pdf,cfd_fres_pdf]);
figure(309); hold off; plot(plt_bins,mri_pdf/max_val,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',line_size);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',line_size);
plot(plt_bins,stb_pdf/max_val,'--','Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
axis([-0.4 0.6 0 1])
xlabel('v velocity (m/s)')
ylabel('PDF')
figure(310); hold off; plot(plt_bins,mri_cdf,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_cdf,'Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_cdf,'Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
xlabel('v velocity (m/s)')
ylabel('CDF')

% For w velocity magnitude
min_bin = min([min(cfd_w_dis),min(stb_w_dis),min(mri_w_dis)]);
max_bin = max([max(cfd_w_dis),max(stb_w_dis),max(mri_w_dis)]);
bin_edges = linspace(min_bin,max_bin,num_bins+1);
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
figure(1);
mri_pdf = histogram(mri_w_dis,bin_edges,'Normalization','pdf');
mri_pdf = mri_pdf.Values;
mri_cdf = histogram(mri_w_dis,bin_edges,'Normalization','cdf');
mri_cdf = mri_cdf.Values;
stb_pdf = histogram(stb_w_dis,bin_edges,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_w_dis,bin_edges,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_fres_w_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_fres_pdf = histogram(cfd_fres_w_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
cfd_pdf = histogram(cfd_w_dis,bin_edges,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_w_dis,bin_edges,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
max_val = max([mri_pdf,stb_pdf,cfd_pdf,stb_fres_pdf,cfd_fres_pdf]);
figure(311); hold off; plot(plt_bins,mri_pdf/max_val,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',line_size);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',line_size);
plot(plt_bins,stb_pdf/max_val,'--','Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
axis([-0.2 0.4 0 1])
xlabel('w velocity (m/s)')
ylabel('PDF')

figure(312); hold off; plot(plt_bins,mri_cdf,'Color',cmri,'LineWidth',line_size);
hold on; plot(plt_bins,stb_cdf,'Color',cstb,'LineWidth',line_size);
plot(plt_bins,cfd_cdf,'Color',ccfd,'LineWidth',line_size);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
xlabel('velocity (m/s)')
ylabel('CDF')

% Print the distribution plots
print(305,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_velmag_pdf.pdf'])
print(306,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_velmag_cdf.pdf'])
print(307,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_uvel_pdf.pdf'])
print(308,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_uvel_cdf.pdf'])
print(309,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_vvel_pdf.pdf'])
print(310,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_vvel_cdf.pdf'])
print(311,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_wvel_pdf.pdf'])
print(312,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_wvel_cdf.pdf'])

% % % figure(303); hold off; histogram(mri_vel_dis,bin_edges,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor',cmri,'LineStyle','none','Normalization','pdf')
% % % hold on; histogram(stb_vel_dis,bin_edges,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor',cstb,'LineStyle','none','Normalization','pdf')
% % % histogram(cfd_vel_dis,bin_edges,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor',ccfd,'LineStyle','none','Normalization','pdf')
% % % mri_pdf = histogram(mri_vel_dis,bin_edges,'Normalization','pdf');
% % % 
% % % figure(304); histogram(mri_vel_dis,bin_edges,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor',cmri,'LineStyle','none','Normalization','cdf')
% % % hold on; histogram(stb_vel_dis,bin_edges,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor',cstb,'LineStyle','none','Normalization','cdf')
% % % histogram(cfd_vel_dis,bin_edges,'EdgeColor','none','EdgeAlpha',0.8,'FaceAlpha',0.8,'FaceColor',ccfd,'LineStyle','none','Normalization','cdf')

w = 1;


end


