function [] = plot_wall_shear_stress_v3(cadat,fslash,PLOT_WSS_SURFACE,PLOT_OSI_SURFACE)
%close all

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
color_order = [cmri;cstb_light;ccfd_light;cstb;ccfd];

% Set file number of peak systole for each modality
if strcmp(cadat.FILE.anynum,'111')
    ps_mri = 14; % For ANY111 = 14, for ANY007 = 4
    ps_stb = 213; % For ANY111 = 213, for ANY007 = 102 
    ps_cfd = 321; % For ANY111 = 323, for ANY007 = 102
    % Set the file numbers
    mri_file_s = 1; mri_file_e = 20;
    stb_file_s = 1; stb_file_e = 304;
    cfd_file_s = 1; cfd_file_e = 541;
else
    ps_mri = 4; % For ANY111 = 14, for ANY007 = 4
    ps_stb = 97; % For ANY111 = 213, for ANY007 = 97
    ps_cfd = 102; % For ANY111 = 323, for ANY007 = 102
    % Set the file numbers
    mri_file_s = 1; mri_file_e = 13;
    stb_file_s = 1; stb_file_e = 365;
    cfd_file_s = 1; cfd_file_e = 386;
end

%%% Load wall shear stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the wall shear stress for peak systole
fp = load([cadat.DIR.saveppfiles,fslash,'wall_shear_stress',fslash,'any',cadat.FILE.anynum,'_MRI_wss_',num2str(ps_mri,'%05i'),'.mat']);
mri_wssmag = fp.wss_mag; % MRI
fp = load([cadat.DIR.saveppfiles,fslash,'wall_shear_stress',fslash,'any',cadat.FILE.anynum,'_STB_wss_',num2str(ps_stb,'%05i'),'.mat']);
stb_wssmag = fp.wss_mag; % Voxel averaged resolution shake the box
fp = load([cadat.DIR.saveppfiles,fslash,'wall_shear_stress',fslash,'any',cadat.FILE.anynum,'_CFD_wss_',num2str(ps_cfd,'%05i'),'.mat']);
cfd_wssmag = fp.wss_mag; % Voxel averaged resolution CFD
fp = load([cadat.DIR.saveppfiles,fslash,'wall_shear_stress',fslash,'any',cadat.FILE.anynum,'_STB_full_resolution_wss_',num2str(ps_stb,'%05i'),'.mat']);
stb_fres_wssmag = fp.wss_mag; % Full resolution shake the box
fp = load([cadat.DIR.saveppfiles,fslash,'wall_shear_stress',fslash,'any',cadat.FILE.anynum,'_CFD_full_resolution_wss_',num2str(ps_cfd,'%05i'),'.mat']);
cfd_fres_wssmag = fp.wss_mag; % Full resolution CFD

%%% Load x,y,z coordinates from velocity mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the MRI velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_MRI_masked_velocity_grid.mat']);
mri_velmask = fp.velmask; % Velocity grid mask
mri_maskX = fp.xmask;
mri_maskY = fp.ymask;
mri_maskZ = fp.zmask;
% Load the STB velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_STB_masked_velocity_grid.mat']);
stb_velmask = fp.velmask; % Velocity grid mask
stb_maskX = fp.xmask;
stb_maskY = fp.ymask;
stb_maskZ = fp.zmask;
% Load the CFD velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_CFD_masked_velocity_grid.mat']);
cfd_velmask = fp.velmask; % Velocity grid mask
cfd_maskX = fp.xmask;
cfd_maskY = fp.ymask;
cfd_maskZ = fp.zmask;
% Load the STB voxel averaged velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_STB_masked_velocity_grid-voxel_averaged.mat']);
stb_velmask_va = fp.velmask; % Velocity grid mask
stb_maskX_va = fp.xmask;
stb_maskY_va = fp.ymask;
stb_maskZ_va = fp.zmask;
% Load the CFD voxel averaged velocity mask
fp = load([cadat.DIR.savefiles,'ANY-',cadat.FILE.anynum,'_CFD_masked_velocity_grid-voxel_averaged.mat']);
cfd_velmask_va = fp.velmask; % Velocity grid mask
cfd_maskX_va = fp.xmask;
cfd_maskY_va = fp.ymask;
cfd_maskZ_va = fp.zmask;

%%% Load the aneurysm masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp = load([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_MRI_aneursym_mask.mat']);
mri_any_mask = fp.any_mask;
fp = load([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_STB_aneursym_mask.mat']);
stb_any_mask = fp.any_mask;
fp = load([cadat.DIR.saveppfiles,'ANY-',cadat.FILE.anynum,'_CFD_aneursym_mask.mat']);
cfd_any_mask = fp.any_mask;

%%% Reshape wall shear stress for plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape wall shear stress vector - MRI
mri_wss = reshape(mri_wssmag,[size(mri_wssmag,1)*size(mri_wssmag,2)*size(mri_wssmag,3),1]);
% Reshape the mask grid points
[xm,ym,zm] = meshgrid(mri_maskX,mri_maskY,mri_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = mri_wss > 0;
mri_wssf = mri_wss(kp_inds==1);
mri_wss_x = xml(kp_inds==1);
mri_wss_y = yml(kp_inds==1);
mri_wss_z = zml(kp_inds==1);

% Reshape wall shear stress vector - STB Full Resolution
stb_fres_wss = reshape(stb_fres_wssmag,[size(stb_fres_wssmag,1)*size(stb_fres_wssmag,2)*size(stb_fres_wssmag,3),1]);
% Reshape the mask grid points
[xm,ym,zm] = meshgrid(stb_maskX,stb_maskY,stb_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = stb_fres_wss > 0;
stb_fres_wssf = stb_fres_wss(kp_inds==1);
stb_fres_wss_x = xml(kp_inds==1);
stb_fres_wss_y = yml(kp_inds==1);
stb_fres_wss_z = zml(kp_inds==1);

% Reshape wall shear stress vector - CFD
cfd_fres_wss = reshape(cfd_fres_wssmag,[size(cfd_fres_wssmag,1)*size(cfd_fres_wssmag,2)*size(cfd_fres_wssmag,3),1]);
% Reshape the mask grid points
[xm,ym,zm] = meshgrid(cfd_maskX,cfd_maskY,cfd_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = cfd_fres_wss > 0;
cfd_fres_wssf = cfd_fres_wss(kp_inds==1);
cfd_fres_wss_x = xml(kp_inds==1);
cfd_fres_wss_y = yml(kp_inds==1);
cfd_fres_wss_z = zml(kp_inds==1);

% Reshape wall shear stress vector - STB voxel averaged
stb_wss = reshape(stb_wssmag,[size(stb_wssmag,1)*size(stb_wssmag,2)*size(stb_wssmag,3),1]);
% Reshape the mask grid points
[xm,ym,zm] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = stb_wss > 0;
stb_wssf = stb_wss(kp_inds==1);
stb_wss_x = xml(kp_inds==1);
stb_wss_y = yml(kp_inds==1);
stb_wss_z = zml(kp_inds==1);

% Reshape wall shear stress vector - CFD voxel averaged
cfd_wss = reshape(cfd_wssmag,[size(cfd_wssmag,1)*size(cfd_wssmag,2)*size(cfd_wssmag,3),1]);
% Reshape the mask grid points
[xm,ym,zm] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = cfd_wss > 0;
cfd_wssf = cfd_wss(kp_inds==1);
cfd_wss_x = xml(kp_inds==1);
cfd_wss_y = yml(kp_inds==1);
cfd_wss_z = zml(kp_inds==1);

%%% Load the STL files for all modalities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vert_stb,face_stb] = stlread([cadat.DIR.stlfiles,'any',cadat.FILE.anynum,'_STB_velocity_plot.stl']);
[vert_cfd,face_cfd] = stlread([cadat.DIR.stlfiles,'any',cadat.FILE.anynum,'_CFD_velocity_plot.stl']);
[vert_mri,face_mri] = stlread([cadat.DIR.stlfiles,'any',cadat.FILE.anynum,'_MRI_velocity_plot.stl']);

%%% Set the axis range
if strcmp(cadat.FILE.anynum,'111')
    %axis_range = [-10 10 -10 15 -5 10];
    %view_angle = [45,50];
    axis_range = [-10 10 -5 10 -5 15];
    view_angle = [-166,8];
    view_angle_cfd = [-166,8];
    max_color = 25;
else
    %axis_range = [-15 15 -8 24 -5 15];
    %view_angle = [-51,48];
    axis_range = [-15 15 -5 15 -8 24];
    view_angle = [158,-3];
    view_angle_cfd = [150,-3];
    max_color = 10;
end

if PLOT_WSS_SURFACE
    %%% Interpolate the color scheme for the WSS contours %%%%%%%%%%%%%%%%%%%%%
    % For each vertex interpolate the WSS value based on the wall shear stress
    % values
    min_pts = 3;

    % Iterate through all modalities and resolutions
    for ftype = 1:1:5
        if ftype == 1
            % Set the current x, y, z, and wss vectors
            curr_x = mri_wss_x;
            curr_y = mri_wss_y;
            curr_z = mri_wss_z;
            curr_wss = mri_wssf;
            curr_vert = vert_mri;
            min_dist = 1.25;
        elseif ftype == 2
            % Set the current x, y, z, and wss vectors
            curr_x = stb_fres_wss_x;
            curr_y = stb_fres_wss_y;
            curr_z = stb_fres_wss_z;
            curr_wss = stb_fres_wssf;
            curr_vert = vert_stb;
            min_dist = 0.45;
        elseif ftype == 3
            % Set the current x, y, z, and wss vectors
            curr_x = cfd_fres_wss_x;
            curr_y = cfd_fres_wss_y;
            curr_z = cfd_fres_wss_z;
            curr_wss = cfd_fres_wssf;
            curr_vert = vert_cfd;
            min_dist = 0.45;
        elseif ftype == 4
            % Set the current x, y, z, and wss vectors
            curr_x = stb_wss_x;
            curr_y = stb_wss_y;
            curr_z = stb_wss_z;
            curr_wss = stb_wssf;
            curr_vert = vert_stb;
            min_dist = 1.25;
        else
            % Set the current x, y, z, and wss vectors
            curr_x = cfd_wss_x;
            curr_y = cfd_wss_y;
            curr_z = cfd_wss_z;
            curr_wss = cfd_wssf;
            curr_vert = vert_cfd;
            min_dist = 1.25;
        end

        % Initialize the output color matrix and interpolate the color
        curr_wss_color = zeros(size(curr_vert,1),1);
        for q = 1:1:size(curr_vert,1)
            % Get the current vertex 
            cx = curr_vert(q,1);
            cy = curr_vert(q,2);
            cz = curr_vert(q,3);

            % Get the closest WSS points to the current vertex
            pt_dists = sqrt((curr_x-cx).^2 + (curr_y-cy).^2 + (curr_z-cz).^2);
            [srt_dists,srt_pts] = sort(pt_dists,'ascend');
            srt_wss = curr_wss(srt_pts);

            % Identify the points which are to be included in the interpolated
            % color calculation
            curr_dist = min_dist; % Set the current distance
            num_pts = sum(srt_dists <= curr_dist); % Compute the number of points in the minimum distance
            while num_pts < min_pts % Iterate until the minimum number of points are included
                curr_dist = curr_dist + 0.1; % Increase the distance 
                num_pts = sum(srt_dists <= curr_dist); % Recompute the number of points included
            end
            % Collect the indices to be included in the calculation
            inc_inds = (srt_dists <= curr_dist); % Gather the indices of points to include
            inc_wss = srt_wss(inc_inds == 1); % Gather the WSS of points to include
            inc_dists = srt_dists(inc_inds == 1); % Gather the distances of points to include
            inc_wgts = 1./(inc_dists.^2); % Compute the weights
            inc_wgts = inc_wgts/sum(inc_wgts); % Normalize the weights
            interp_wss = sum(inc_wgts.*inc_wss);

            % Save the interpolated WSS in the wss color vector
            curr_wss_color(q) = interp_wss;
        end
        % Set the current color to the appropriate array
        if ftype == 1
            mri_wss_color = curr_wss_color;
        elseif ftype == 2
            stb_fres_wss_color = curr_wss_color; 
        elseif ftype == 3
            cfd_fres_wss_color = curr_wss_color;
        elseif ftype == 4
            stb_wss_color = curr_wss_color; 
        else
            cfd_wss_color = curr_wss_color;
        end
    end

    %%% Plot the WSS contours %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = load('wss_cmap.mat');
    wss_cmap = fp.cmap;
    fsize = 28;

    % Plot the surface contours
    figure(440); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_mri,'Faces',face_mri,'FaceVertexCData',10*mri_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('z (mm)'); zlabel('y (mm)'); colorbar;

    figure(441); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_stb,'Faces',face_stb,'FaceVertexCData',10*stb_fres_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('z (mm)'); zlabel('y (mm)'); colorbar;

    figure(442); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_cfd,'Faces',face_cfd,'FaceVertexCData',10*cfd_fres_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;

    figure(443); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_stb,'Faces',face_stb,'FaceVertexCData',10*stb_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;

    figure(444); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_cfd,'Faces',face_cfd,'FaceVertexCData',10*cfd_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;


    %print(440,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_mri_peaksystole_wss_contour.pdf'])
    %print(441,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_fres_peaksystole_wss_contour.pdf'])
    %print(442,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_fres_peaksystole_wss_contour.pdf'])
    %print(443,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_peaksystole_wss_contour.pdf'])
    %print(444,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_peaksystole_wss_contour.pdf'])


    %%% Plot the normalized WSS contours %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the max to ensure a single max is not used
    srt_wss = sort(mri_wss_color,'descend');
    max_mri = median(srt_wss(1:5));
    srt_wss = sort(stb_fres_wss_color,'descend');
    max_stb_fres = median(srt_wss(1:5));
    srt_wss = sort(cfd_fres_wss_color,'descend');
    max_cfd_fres = median(srt_wss(1:5));
    srt_wss = sort(stb_wss_color,'descend');
    max_stb = median(srt_wss(1:5));
    srt_wss = sort(cfd_wss_color,'descend');
    max_cfd = median(srt_wss(1:5));

    fsize = 28;
    max_color = max_color/100;
    
    % Flip y and z for better plotting
    vert_mri2 = [vert_mri(:,1),vert_mri(:,3),vert_mri(:,2)];
    face_mri2 = [face_mri(:,1),face_mri(:,3),face_mri(:,2)];
    vert_stb2 = [vert_stb(:,1),vert_stb(:,3),vert_stb(:,2)];
    face_stb2 = [face_stb(:,1),face_stb(:,3),face_stb(:,2)];
    vert_cfd2 = [vert_cfd(:,1),vert_cfd(:,3),vert_cfd(:,2)];
    face_cfd2 = [face_cfd(:,1),face_cfd(:,3),face_cfd(:,2)];
    
    % Plot the surface contours
    figure(430); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_mri2,'Faces',face_mri2,'FaceVertexCData',mri_wss_color/max_mri,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); grid(gca,'on')
    set(gca,'XDir','reverse');
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('z (mm)'); zlabel('y (mm)'); colorbar;
    axis off;

    figure(431); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_stb2,'Faces',face_stb2,'FaceVertexCData',stb_fres_wss_color/max_stb_fres,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); grid(gca,'on')
    set(gca,'XDir','reverse');
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('z (mm)'); zlabel('y (mm)'); colorbar;
    axis off;

    figure(432); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_cfd2,'Faces',face_cfd2,'FaceVertexCData',cfd_fres_wss_color/max_cfd_fres,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); grid(gca,'on')
    set(gca,'XDir','reverse');
    axis(axis_range); view(view_angle_cfd); xlabel('x (mm)'); ylabel('z (mm)'); zlabel('y (mm)'); colorbar;
    axis off;

    figure(433); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_stb2,'Faces',face_stb2,'FaceVertexCData',stb_wss_color/max_stb,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); grid(gca,'on')
    set(gca,'XDir','reverse');
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('z (mm)'); zlabel('y (mm)'); colorbar;
    axis off;

    figure(434); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_cfd2,'Faces',face_cfd2,'FaceVertexCData',cfd_wss_color/max_cfd,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); grid(gca,'on')
    set(gca,'XDir','reverse');
    axis(axis_range); view(view_angle_cfd); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); h = colorbar;
    axis off;

    print(430,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_normalized_mri_peaksystole_wss_contour.pdf'])
    print(431,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_normalized_stb_fres_peaksystole_wss_contour.pdf'])
    print(432,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_normalized_cfd_fres_peaksystole_wss_contour.pdf'])
    print(433,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_normalized_stb_peaksystole_wss_contour.pdf'])
    print(434,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_normalized_cfd_peaksystole_wss_contour.pdf'])
    
    % ylabel(h,'normalized wall shear stress');
    % print(434,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_normalized_wss_contour_colorbar.pdf'])
end

%%% Plot the scatter points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure(449); hold off; scatter3(mri_wss_x,mri_wss_y,mri_wss_z,10,10*mri_wssf,'filled')
% % caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); colormap(gca,wss_cmap);
% % axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
% % figure(450); hold off; scatter3(stb_wss_x,stb_wss_y,stb_wss_z,10,10*stb_wssf,'filled')
% % caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); colormap(gca,wss_cmap);
% % axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
% % figure(451); hold off; scatter3(cfd_wss_x,cfd_wss_y,cfd_wss_z,10,10*cfd_wssf,'filled');
% % caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); colormap(gca,wss_cmap);
% % axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
% % figure(452); hold off; scatter3(stb_fres_wss_x,stb_fres_wss_y,stb_fres_wss_z,10,10*stb_fres_wssf,'filled');
% % caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); colormap(gca,wss_cmap);
% % axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
% % figure(453); hold off; scatter3(cfd_fres_wss_x,cfd_fres_wss_y,cfd_fres_wss_z,10,10*cfd_fres_wssf,'filled');
% % caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); colormap(gca,wss_cmap);
% % axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
% % print(449,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_mri_peaksystole_wss.pdf'])
% % print(450,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_peaksystole_wss.pdf'])
% % print(451,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_peaksystole_wss.pdf'])
% % print(452,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_fres_peaksystole_wss.pdf'])
% % print(453,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_fres_peaksystole_wss.pdf'])


%%% Get wall shear stress statistics in aneurysm %%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the wall shear stress in the aneurysm
% MRI mask velocity points
any_inds = zeros(size(mri_wss_x));
% Load STB velocity field
for q = 1:1:size(mri_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = mri_any_mask(q,1);
    cy = mri_any_mask(q,2);
    cz = mri_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((mri_wss_x-cx).^2 + (mri_wss_y-cy).^2 + (mri_wss_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
mri_wss_dis = mri_wssf(any_inds==1); 

% STB mask velocity points (STB voxel averaged resolution)
any_inds = zeros(size(stb_wss_x));
% Load STB velocity field
for q = 1:1:size(stb_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = stb_any_mask(q,1);
    cy = stb_any_mask(q,2);
    cz = stb_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((stb_wss_x-cx).^2 + (stb_wss_y-cy).^2 + (stb_wss_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
stb_wss_dis = stb_wssf(any_inds==1); 

% CFD mask velocity points (CFD voxel averaged resolution)
any_inds = zeros(size(cfd_wss_x));
% Load STB velocity field
for q = 1:1:size(cfd_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = cfd_any_mask(q,1);
    cy = cfd_any_mask(q,2);
    cz = cfd_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((cfd_wss_x-cx).^2 + (cfd_wss_y-cy).^2 + (cfd_wss_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
cfd_wss_dis = cfd_wssf(any_inds==1);

% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
stb_vamask = reshape(stb_velmask_va,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskX = reshape(Xva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskY = reshape(Yva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskZ = reshape(Zva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(stb_fres_wss_x));
for q = 1:1:size(stb_fres_wss_x)
    % Get current x,y,z point
    cx = stb_fres_wss_x(q);
    cy = stb_fres_wss_y(q);
    cz = stb_fres_wss_z(q);
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
stb_wss_fres_dis = stb_fres_wssf(fres_inany==1);


% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
cfd_vamask = reshape(cfd_velmask_va,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskX = reshape(Xva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskY = reshape(Yva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskZ = reshape(Zva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(cfd_fres_wss_x));
for q = 1:1:size(cfd_fres_wss_x)
    % Get current x,y,z point
    cx = cfd_fres_wss_x(q);
    cy = cfd_fres_wss_y(q);
    cz = cfd_fres_wss_z(q);
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
cfd_wss_fres_dis = cfd_fres_wssf(fres_inany==1);


% Remove any wall shear stress points whose value is greater than 2
% standard deviations away
stb_rm = stb_wss_dis < (mean(stb_wss_dis)+2*std(stb_wss_dis));
stb_wss_dis = 10*stb_wss_dis;%(stb_rm); % Convert to dynes/cm^2
cfd_rm = cfd_wss_dis < (mean(cfd_wss_dis)+2*std(cfd_wss_dis));
cfd_wss_dis = 10*cfd_wss_dis;%(cfd_rm); % Convert to dynes/cm^2
stb_wss_fres_dis = 10*stb_wss_fres_dis;
cfd_wss_fres_dis = 10*cfd_wss_fres_dis;
mri_wss_dis = 10*mri_wss_dis;

% Plot the histogram of points
% Voxel average bins
num_bins = round(length(mri_wss_dis)/10);%8;
min_bin_va = min([min(cfd_wss_dis),min(stb_wss_dis),min(mri_wss_dis),min(stb_wss_fres_dis),min(cfd_wss_fres_dis)]);
max_bin_va = max([max(cfd_wss_dis),max(stb_wss_dis),max(mri_wss_dis),max(stb_wss_fres_dis),max(cfd_wss_fres_dis)]);
last_bin_edge = 15;%ceil(max_bin_va);
bin_edges_va = linspace(0,last_bin_edge,num_bins+1);
%bin_edges_va = [bin_edges_va,max_bin_va];
plt_bins_va = (bin_edges_va(2:end)-bin_edges_va(1:end-1))/2 + bin_edges_va(1:end-1);
% Full resolution bins
num_bins = round(length(cfd_wss_fres_dis)/50);%25;
min_bin = min([min(cfd_wss_dis),min(stb_wss_dis),min(stb_wss_fres_dis),min(cfd_wss_fres_dis),min(mri_wss_dis)]);
max_bin = max([max(cfd_wss_dis),max(stb_wss_dis),max(stb_wss_fres_dis),max(cfd_wss_fres_dis),max(mri_wss_dis)]);
last_bin_edge = 15;%ceil(max_bin);
bin_edges = linspace(0,last_bin_edge,num_bins+1);
%bin_edges = [bin_edges,max_bin];
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
% Compute histogram
figure(1);
stb_pdf = histogram(stb_wss_dis,bin_edges_va,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_wss_dis,bin_edges_va,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_wss_fres_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_pdf = histogram(cfd_wss_dis,bin_edges_va,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_wss_dis,bin_edges_va,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
cfd_fres_pdf = histogram(cfd_wss_fres_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
mri_pdf = histogram(mri_wss_dis,bin_edges_va,'Normalization','pdf');
mri_pdf = mri_pdf.Values;


% Plot the PDFS
lsize = 3;
fsize = 24;
max_val = max([mri_pdf,cfd_pdf,stb_pdf,cfd_fres_pdf,stb_fres_pdf]);
figure(511); hold off; plot(plt_bins_va,mri_pdf/max_val,'Color',cmri,'LineWidth',lsize);
hold on; plot(plt_bins_va,stb_pdf/max_val,'--','Color',cstb,'LineWidth',lsize);
plot(plt_bins_va,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',lsize);
plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',lsize);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',lsize);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); set(gca,'Box','off')
xlabel('wall shear stress (dynes/cm^{2})')
ylabel('PDF')
pbaspect([1 0.6 1])
% figure(512); hold off; plot(plt_bins_va,stb_cdf,'Color',cstb,'LineWidth',1.5);
% hold on; plot(plt_bins_va,cfd_cdf,'Color',ccfd,'LineWidth',1.5);
% set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);  set(gca,'Box','off')
% xlabel('wall shear stress (Pa)')
% ylabel('CDF')

print(511,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_wss_pdf.pdf'])
savefig(511,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_wss_pdf.fig'])
%print(512,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_wss_cdf.pdf'])


%%% Time averaged wall shear stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load time averaged wall shear stress for each case
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_MRI_TAWSS','.mat']);
mri_tawss = fp.tawss_mag;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_STB_TAWSS','.mat']);
stb_tawss = fp.tawss_mag;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_CFD_TAWSS','.mat']);
cfd_tawss = fp.tawss_mag;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_STB_full_resolution_TAWSS','.mat']);
stb_fres_tawss = fp.tawss_mag;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_CFD_full_resolution_TAWSS','.mat']);
cfd_fres_tawss = fp.tawss_mag;

% Reshape the TAWSS to line arrays
mri_tawss = reshape(mri_tawss,[size(mri_tawss,1)*size(mri_tawss,2)*size(mri_tawss,3),1]);
[xm,ym,zm] = meshgrid(mri_maskX,mri_maskY,mri_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = mri_tawss > 0;
mri_tawss = mri_tawss(kp_inds==1);
mri_tawss_x = xml(kp_inds==1);
mri_tawss_y = yml(kp_inds==1);
mri_tawss_z = zml(kp_inds==1);
stb_tawss = reshape(stb_tawss,[size(stb_tawss,1)*size(stb_tawss,2)*size(stb_tawss,3),1]);
[xm,ym,zm] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = stb_tawss > 0;
stb_tawss = stb_tawss(kp_inds==1);
stb_tawss_x = xml(kp_inds==1);
stb_tawss_y = yml(kp_inds==1);
stb_tawss_z = zml(kp_inds==1);
cfd_tawss = reshape(cfd_tawss,[size(cfd_tawss,1)*size(cfd_tawss,2)*size(cfd_tawss,3),1]);
[xm,ym,zm] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = cfd_tawss > 0;
cfd_tawss = cfd_tawss(kp_inds==1);
cfd_tawss_x = xml(kp_inds==1);
cfd_tawss_y = yml(kp_inds==1);
cfd_tawss_z = zml(kp_inds==1);
stb_fres_tawss = reshape(stb_fres_tawss,[size(stb_fres_tawss,1)*size(stb_fres_tawss,2)*size(stb_fres_tawss,3),1]);
[xm,ym,zm] = meshgrid(stb_maskX,stb_maskY,stb_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = stb_fres_tawss > 0;
stb_fres_tawss = stb_fres_tawss(kp_inds==1);
stb_fres_tawss_x = xml(kp_inds==1);
stb_fres_tawss_y = yml(kp_inds==1);
stb_fres_tawss_z = zml(kp_inds==1);
cfd_fres_tawss = reshape(cfd_fres_tawss,[size(cfd_fres_tawss,1)*size(cfd_fres_tawss,2)*size(cfd_fres_tawss,3),1]);
[xm,ym,zm] = meshgrid(cfd_maskX,cfd_maskY,cfd_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = cfd_fres_tawss > 0;
cfd_fres_tawss = cfd_fres_tawss(kp_inds==1);
cfd_fres_tawss_x = xml(kp_inds==1);
cfd_fres_tawss_y = yml(kp_inds==1);
cfd_fres_tawss_z = zml(kp_inds==1);

% Save the TAWSS values into an array for Bland-Altman analysis
tawss_all{1,1} = mri_tawss;
tawss_all_X{1,1} = mri_tawss_x;
tawss_all_Y{1,1} = mri_tawss_y;
tawss_all_Z{1,1} = mri_tawss_z;
tawss_all{2,1} = stb_tawss;
tawss_all_X{2,1} = stb_tawss_x;
tawss_all_Y{2,1} = stb_tawss_y;
tawss_all_Z{2,1} = stb_tawss_z;
tawss_all{3,1} = cfd_tawss;
tawss_all_X{3,1} = cfd_tawss_x;
tawss_all_Y{3,1} = cfd_tawss_y;
tawss_all_Z{3,1} = cfd_tawss_z;
tawss_all{4,1} = stb_fres_tawss;
tawss_all_X{4,1} = stb_fres_tawss_x;
tawss_all_Y{4,1} = stb_fres_tawss_y;
tawss_all_Z{4,1} = stb_fres_tawss_z;
tawss_all{5,1} = cfd_fres_tawss;
tawss_all_X{5,1} = cfd_fres_tawss_x;
tawss_all_Y{5,1} = cfd_fres_tawss_y;
tawss_all_Z{5,1} = cfd_fres_tawss_z;

%%% Get TAWSS in the aneurysm only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the wall shear stress in the aneurysm
% MRI mask velocity points
any_inds = zeros(size(mri_tawss_x));
% Load STB velocity field
for q = 1:1:size(mri_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = mri_any_mask(q,1);
    cy = mri_any_mask(q,2);
    cz = mri_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((mri_tawss_x-cx).^2 + (mri_tawss_y-cy).^2 + (mri_tawss_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
mri_tawss_dis = mri_tawss(any_inds==1); 

% STB mask velocity points (STB voxel averaged resolution)
any_inds = zeros(size(stb_tawss_x));
% Load STB velocity field
for q = 1:1:size(stb_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = stb_any_mask(q,1);
    cy = stb_any_mask(q,2);
    cz = stb_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((stb_tawss_x-cx).^2 + (stb_tawss_y-cy).^2 + (stb_tawss_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
stb_tawss_dis = stb_tawss(any_inds==1); 

% CFD mask velocity points (CFD voxel averaged resolution)
any_inds = zeros(size(cfd_tawss_x));
% Load STB velocity field
for q = 1:1:size(cfd_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = cfd_any_mask(q,1);
    cy = cfd_any_mask(q,2);
    cz = cfd_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((cfd_tawss_x-cx).^2 + (cfd_tawss_y-cy).^2 + (cfd_tawss_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
cfd_tawss_dis = cfd_tawss(any_inds==1);

% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
stb_vamask = reshape(stb_velmask_va,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskX = reshape(Xva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskY = reshape(Yva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskZ = reshape(Zva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(stb_fres_tawss_x));
for q = 1:1:size(stb_fres_tawss_x)
    % Get current x,y,z point
    cx = stb_fres_tawss_x(q);
    cy = stb_fres_tawss_y(q);
    cz = stb_fres_tawss_z(q);
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
stb_tawss_fres_dis = stb_fres_tawss(fres_inany==1);


% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
cfd_vamask = reshape(cfd_velmask_va,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskX = reshape(Xva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskY = reshape(Yva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskZ = reshape(Zva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(cfd_fres_tawss_x));
for q = 1:1:size(cfd_fres_tawss_x)
    % Get current x,y,z point
    cx = cfd_fres_tawss_x(q);
    cy = cfd_fres_tawss_y(q);
    cz = cfd_fres_tawss_z(q);
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
cfd_tawss_fres_dis = cfd_fres_tawss(fres_inany==1);


% Plot the TAWSS
fp = load('wss_cmap.mat');
wss_cmap = fp.cmap;
max_color = 10;
fsize = 20;
figure(459); hold off; scatter3(mri_tawss_x,mri_tawss_y,mri_tawss_z,10,10*mri_tawss,'filled')
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(460); hold off; scatter3(stb_tawss_x,stb_tawss_y,stb_tawss_z,10,10*stb_tawss,'filled')
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(461); hold off; scatter3(cfd_tawss_x,cfd_tawss_y,cfd_tawss_z,10,10*cfd_tawss,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(462); hold off; scatter3(stb_fres_tawss_x,stb_fres_tawss_y,stb_fres_tawss_z,10,10*stb_fres_tawss,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(463); hold off; scatter3(cfd_fres_tawss_x,cfd_fres_tawss_y,cfd_fres_tawss_z,10,10*cfd_fres_tawss,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
print(459,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_mri_tawss.pdf'])
print(460,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_tawss.pdf'])
print(461,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_tawss.pdf'])
print(462,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_fres_tawss.pdf'])
print(463,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_fres_tawss.pdf'])

% Compute time and space averaged TAWSS
mri_tawss_avg = mean(mri_tawss_dis);
stb_tawss_avg = mean(stb_tawss_dis);
cfd_tawss_avg = mean(cfd_tawss_dis);
stb_fres_tawss_avg = mean(stb_tawss_fres_dis);
cfd_fres_tawss_avg = mean(cfd_tawss_fres_dis);
tawss_avg = 10*[mri_tawss_avg,stb_fres_tawss_avg,cfd_fres_tawss_avg,stb_tawss_avg,cfd_tawss_avg];

%%% Time averaged OSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load OSI for each case
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_MRI_OSI','.mat']);
mri_osi = fp.osi_avg;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_STB_OSI','.mat']);
stb_osi = fp.osi_avg;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_CFD_OSI','.mat']);
cfd_osi = fp.osi_avg;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_STB_full_resolution_OSI','.mat']);
stb_fres_osi = fp.osi_avg;
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_CFD_full_resolution_OSI','.mat']);
cfd_fres_osi = fp.osi_avg;

% Reshape the OSI to line arrays
mri_osi = reshape(mri_osi,[size(mri_osi,1)*size(mri_osi,2)*size(mri_osi,3),1]);
[xm,ym,zm] = meshgrid(mri_maskX,mri_maskY,mri_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = mri_osi > 0;
mri_osi = mri_osi(kp_inds==1);
mri_osi_x = xml(kp_inds==1);
mri_osi_y = yml(kp_inds==1);
mri_osi_z = zml(kp_inds==1);
stb_osi = reshape(stb_osi,[size(stb_osi,1)*size(stb_osi,2)*size(stb_osi,3),1]);
[xm,ym,zm] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = stb_osi >= 0;
stb_osi = stb_osi(kp_inds==1);
stb_osi_x = xml(kp_inds==1);
stb_osi_y = yml(kp_inds==1);
stb_osi_z = zml(kp_inds==1);
cfd_osi = reshape(cfd_osi,[size(cfd_osi,1)*size(cfd_osi,2)*size(cfd_osi,3),1]);
[xm,ym,zm] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = cfd_osi >= 0;
cfd_osi = cfd_osi(kp_inds==1);
cfd_osi_x = xml(kp_inds==1);
cfd_osi_y = yml(kp_inds==1);
cfd_osi_z = zml(kp_inds==1);
stb_fres_osi = reshape(stb_fres_osi,[size(stb_fres_osi,1)*size(stb_fres_osi,2)*size(stb_fres_osi,3),1]);
[xm,ym,zm] = meshgrid(stb_maskX,stb_maskY,stb_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = stb_fres_osi >= 0;
stb_fres_osi = stb_fres_osi(kp_inds==1);
stb_fres_osi_x = xml(kp_inds==1);
stb_fres_osi_y = yml(kp_inds==1);
stb_fres_osi_z = zml(kp_inds==1);
cfd_fres_osi = reshape(cfd_fres_osi,[size(cfd_fres_osi,1)*size(cfd_fres_osi,2)*size(cfd_fres_osi,3),1]);
[xm,ym,zm] = meshgrid(cfd_maskX,cfd_maskY,cfd_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = cfd_fres_osi >= 0;
cfd_fres_osi = cfd_fres_osi(kp_inds==1);
cfd_fres_osi_x = xml(kp_inds==1);
cfd_fres_osi_y = yml(kp_inds==1);
cfd_fres_osi_z = zml(kp_inds==1);

% Save the RRT values into an array for Bland-Altman analysis
osi_all{1,1} = mri_osi;
osi_all_X{1,1} = mri_osi_x;
osi_all_Y{1,1} = mri_osi_y;
osi_all_Z{1,1} = mri_osi_z;
osi_all{2,1} = stb_osi;
osi_all_X{2,1} = stb_osi_x;
osi_all_Y{2,1} = stb_osi_y;
osi_all_Z{2,1} = stb_osi_z;
osi_all{3,1} = cfd_osi;
osi_all_X{3,1} = cfd_osi_x;
osi_all_Y{3,1} = cfd_osi_y;
osi_all_Z{3,1} = cfd_osi_z;
osi_all{4,1} = stb_fres_osi;
osi_all_X{4,1} = stb_fres_osi_x;
osi_all_Y{4,1} = stb_fres_osi_y;
osi_all_Z{4,1} = stb_fres_osi_z;
osi_all{5,1} = cfd_fres_osi;
osi_all_X{5,1} = cfd_fres_osi_x;
osi_all_Y{5,1} = cfd_fres_osi_y;
osi_all_Z{5,1} = cfd_fres_osi_z;

%%% Get TAWSS in the aneurysm only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the wall shear stress in the aneurysm
% MRI mask velocity points
any_inds = zeros(size(mri_osi_x));
% Load STB velocity field
for q = 1:1:size(mri_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = mri_any_mask(q,1);
    cy = mri_any_mask(q,2);
    cz = mri_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((mri_osi_x-cx).^2 + (mri_osi_y-cy).^2 + (mri_osi_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
mri_osi_dis = mri_osi(any_inds==1); 

% STB mask velocity points (STB voxel averaged resolution)
any_inds = zeros(size(stb_osi_x));
% Load STB velocity field
for q = 1:1:size(stb_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = stb_any_mask(q,1);
    cy = stb_any_mask(q,2);
    cz = stb_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((stb_osi_x-cx).^2 + (stb_osi_y-cy).^2 + (stb_osi_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
stb_osi_dis = stb_osi(any_inds==1); 

% CFD mask velocity points (CFD voxel averaged resolution)
any_inds = zeros(size(cfd_osi_x));
% Load STB velocity field
for q = 1:1:size(cfd_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = cfd_any_mask(q,1);
    cy = cfd_any_mask(q,2);
    cz = cfd_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((cfd_osi_x-cx).^2 + (cfd_osi_y-cy).^2 + (cfd_osi_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
cfd_osi_dis = cfd_osi(any_inds==1);

% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
stb_vamask = reshape(stb_velmask_va,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskX = reshape(Xva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskY = reshape(Yva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskZ = reshape(Zva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(stb_fres_osi_x));
for q = 1:1:size(stb_fres_osi_x)
    % Get current x,y,z point
    cx = stb_fres_osi_x(q);
    cy = stb_fres_osi_y(q);
    cz = stb_fres_osi_z(q);
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
stb_osi_fres_dis = stb_fres_osi(fres_inany==1);


% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
cfd_vamask = reshape(cfd_velmask_va,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskX = reshape(Xva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskY = reshape(Yva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskZ = reshape(Zva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(cfd_fres_osi_x));
for q = 1:1:size(cfd_fres_osi_x)
    % Get current x,y,z point
    cx = cfd_fres_osi_x(q);
    cy = cfd_fres_osi_y(q);
    cz = cfd_fres_osi_z(q);
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
cfd_osi_fres_dis = cfd_fres_osi(fres_inany==1);

% Plot the osi
fp = load('wss_cmap.mat');
wss_cmap = fp.cmap;
max_color = 0.5;
fsize = 20;
figure(469); hold off; scatter3(mri_osi_x,mri_osi_y,mri_osi_z,10,mri_osi,'filled')
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); %colormap(gca,osi_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(470); hold off; scatter3(stb_osi_x,stb_osi_y,stb_osi_z,10,stb_osi,'filled')
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); %colormap(gca,osi_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(471); hold off; scatter3(cfd_osi_x,cfd_osi_y,cfd_osi_z,10,cfd_osi,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); %colormap(gca,osi_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(472); hold off; scatter3(stb_fres_osi_x,stb_fres_osi_y,stb_fres_osi_z,10,stb_fres_osi,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); %colormap(gca,osi_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(473); hold off; scatter3(cfd_fres_osi_x,cfd_fres_osi_y,cfd_fres_osi_z,10,cfd_fres_osi,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); %colormap(gca,osi_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
print(469,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_mri_osi.pdf'])
print(470,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_osi.pdf'])
print(471,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_osi.pdf'])
print(472,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_fres_osi.pdf'])
print(473,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_fres_osi.pdf'])

% Compute time and space averaged TAWSS
mri_osi_avg = mean(mri_osi_dis);
stb_osi_avg = mean(stb_osi_dis);
cfd_osi_avg = mean(cfd_osi_dis);
stb_fres_osi_avg = mean(stb_osi_fres_dis);
cfd_fres_osi_avg = mean(cfd_osi_fres_dis);
mri_osi_std = std(mri_osi_dis);
stb_osi_std = std(stb_osi_dis);
cfd_osi_std = std(cfd_osi_dis);
stb_fres_osi_std = std(stb_osi_fres_dis);
cfd_fres_osi_std = std(cfd_osi_fres_dis);
osi_avg = [mri_osi_avg,stb_fres_osi_avg,cfd_fres_osi_avg,stb_osi_avg,cfd_osi_avg];
osi_std = [mri_osi_std,stb_fres_osi_std,cfd_fres_osi_std,stb_osi_std,cfd_osi_std];

%%% Time averaged RRT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load RRT for each case
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_MRI_RRT','.mat']);
mri_rrt = fp.rrt_avg; %sqrt(fp.rrt_x.^2 + fp.rrt_y.^2 + fp.rrt_z.^2);
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_STB_RRT','.mat']);
stb_rrt = fp.rrt_avg; %sqrt(fp.rrt_x.^2 + fp.rrt_y.^2 + fp.rrt_z.^2);
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_CFD_RRT','.mat']);
cfd_rrt = fp.rrt_avg; %sqrt(fp.rrt_x.^2 + fp.rrt_y.^2 + fp.rrt_z.^2);
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_STB_full_resolution_RRT','.mat']);
stb_fres_rrt = fp.rrt_avg; %sqrt(fp.rrt_x.^2 + fp.rrt_y.^2 + fp.rrt_z.^2);
fp = load([cadat.DIR.saveppfiles,fslash,'any',cadat.FILE.anynum,'_CFD_full_resolution_RRT','.mat']);
cfd_fres_rrt = fp.rrt_avg; %sqrt(fp.rrt_x.^2 + fp.rrt_y.^2 + fp.rrt_z.^2);

% Reshape the RRT to line arrays
mri_rrt = reshape(mri_rrt,[size(mri_rrt,1)*size(mri_rrt,2)*size(mri_rrt,3),1]);
[xm,ym,zm] = meshgrid(mri_maskX,mri_maskY,mri_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = ~isnan(mri_rrt) & (mri_rrt < 500);
mri_rrt = mri_rrt(kp_inds==1);
mri_rrt_x = xml(kp_inds==1);
mri_rrt_y = yml(kp_inds==1);
mri_rrt_z = zml(kp_inds==1);
stb_rrt = reshape(stb_rrt,[size(stb_rrt,1)*size(stb_rrt,2)*size(stb_rrt,3),1]);
[xm,ym,zm] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = ~isnan(stb_rrt) & (stb_rrt < 500);
stb_rrt = stb_rrt(kp_inds==1);
stb_rrt_x = xml(kp_inds==1);
stb_rrt_y = yml(kp_inds==1);
stb_rrt_z = zml(kp_inds==1);
cfd_rrt = reshape(cfd_rrt,[size(cfd_rrt,1)*size(cfd_rrt,2)*size(cfd_rrt,3),1]);
[xm,ym,zm] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = ~isnan(cfd_rrt) & (cfd_rrt < 500);
cfd_rrt = cfd_rrt(kp_inds==1);
cfd_rrt_x = xml(kp_inds==1);
cfd_rrt_y = yml(kp_inds==1);
cfd_rrt_z = zml(kp_inds==1);
stb_fres_rrt = reshape(stb_fres_rrt,[size(stb_fres_rrt,1)*size(stb_fres_rrt,2)*size(stb_fres_rrt,3),1]);
[xm,ym,zm] = meshgrid(stb_maskX,stb_maskY,stb_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = ~isnan(stb_fres_rrt) & (stb_fres_rrt < 500);
stb_fres_rrt = stb_fres_rrt(kp_inds==1);
stb_fres_rrt_x = xml(kp_inds==1);
stb_fres_rrt_y = yml(kp_inds==1);
stb_fres_rrt_z = zml(kp_inds==1);
cfd_fres_rrt = reshape(cfd_fres_rrt,[size(cfd_fres_rrt,1)*size(cfd_fres_rrt,2)*size(cfd_fres_rrt,3),1]);
[xm,ym,zm] = meshgrid(cfd_maskX,cfd_maskY,cfd_maskZ);
xml = reshape(xm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
yml = reshape(ym,[size(xm,1)*size(xm,2)*size(xm,3),1]);
zml = reshape(zm,[size(xm,1)*size(xm,2)*size(xm,3),1]);
kp_inds = ~isnan(cfd_fres_rrt) & (cfd_fres_rrt < 500);
cfd_fres_rrt = cfd_fres_rrt(kp_inds==1);
cfd_fres_rrt_x = xml(kp_inds==1);
cfd_fres_rrt_y = yml(kp_inds==1);
cfd_fres_rrt_z = zml(kp_inds==1);

% Find the elbow of the RRT values and discard the "exploding" values
% MRI
[srt_rrt,sinds] = sort(mri_rrt);
d2rrt = abs(diff(diff(srt_rrt)));
ct_ind = find(d2rrt > mean(d2rrt),1,'first');
kp_inds = sort(sinds(1:ct_ind),'ascend');
mri_rrt = mri_rrt(kp_inds);
mri_rrt_x = mri_rrt_x(kp_inds);
mri_rrt_y = mri_rrt_y(kp_inds);
mri_rrt_z = mri_rrt_z(kp_inds);
% STB voxel averaged resolution
[srt_rrt,sinds] = sort(stb_rrt);
d2rrt = abs(diff(diff(srt_rrt)));
ct_ind = find(d2rrt > mean(d2rrt),1,'first');
kp_inds = sort(sinds(1:ct_ind),'ascend');
stb_rrt = stb_fres_rrt(kp_inds);
stb_rrt_x = stb_rrt_x(kp_inds);
stb_rrt_y = stb_rrt_y(kp_inds);
stb_rrt_z = stb_rrt_z(kp_inds);
% STB full resolution
[srt_rrt,sinds] = sort(stb_fres_rrt);
d2rrt = abs(diff(diff(srt_rrt)));
ct_ind = find(d2rrt > mean(d2rrt),1,'first');
kp_inds = sort(sinds(1:ct_ind),'ascend');
stb_fres_rrt = stb_fres_rrt(kp_inds);
stb_fres_rrt_x = stb_fres_rrt_x(kp_inds);
stb_fres_rrt_y = stb_fres_rrt_y(kp_inds);
stb_fres_rrt_z = stb_fres_rrt_z(kp_inds);
% CFD voxel averaged resolution
[srt_rrt,sinds] = sort(cfd_rrt);
d2rrt = abs(diff(diff(srt_rrt)));
ct_ind = find(d2rrt > mean(d2rrt),1,'first');
kp_inds = sort(sinds(1:ct_ind),'ascend');
cfd_rrt = cfd_fres_rrt(kp_inds);
cfd_rrt_x = cfd_rrt_x(kp_inds);
cfd_rrt_y = cfd_rrt_y(kp_inds);
cfd_rrt_z = cfd_rrt_z(kp_inds);
% CFD full resolution
[srt_rrt,sinds] = sort(cfd_fres_rrt);
d2rrt = abs(diff(diff(srt_rrt)));
ct_ind = find(d2rrt > mean(d2rrt),1,'first');
kp_inds = sort(sinds(1:ct_ind),'ascend');
cfd_fres_rrt = cfd_fres_rrt(kp_inds);
cfd_fres_rrt_x = cfd_fres_rrt_x(kp_inds);
cfd_fres_rrt_y = cfd_fres_rrt_y(kp_inds);
cfd_fres_rrt_z = cfd_fres_rrt_z(kp_inds);

% Save the RRT values into an array for Bland-Altman analysis
rrt_all{1,1} = mri_rrt;
rrt_all_X{1,1} = mri_rrt_x;
rrt_all_Y{1,1} = mri_rrt_y;
rrt_all_Z{1,1} = mri_rrt_z;
rrt_all{2,1} = stb_rrt;
rrt_all_X{2,1} = stb_rrt_x;
rrt_all_Y{2,1} = stb_rrt_y;
rrt_all_Z{2,1} = stb_rrt_z;
rrt_all{3,1} = cfd_rrt;
rrt_all_X{3,1} = cfd_rrt_x;
rrt_all_Y{3,1} = cfd_rrt_y;
rrt_all_Z{3,1} = cfd_rrt_z;
rrt_all{4,1} = stb_fres_rrt;
rrt_all_X{4,1} = stb_fres_rrt_x;
rrt_all_Y{4,1} = stb_fres_rrt_y;
rrt_all_Z{4,1} = stb_fres_rrt_z;
rrt_all{5,1} = cfd_fres_rrt;
rrt_all_X{5,1} = cfd_fres_rrt_x;
rrt_all_Y{5,1} = cfd_fres_rrt_y;
rrt_all_Z{5,1} = cfd_fres_rrt_z;

%%% Get TAWSS in the aneurysm only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the wall shear stress in the aneurysm
% MRI mask velocity points
any_inds = zeros(size(mri_rrt_x));
% Load STB velocity field
for q = 1:1:size(mri_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = mri_any_mask(q,1);
    cy = mri_any_mask(q,2);
    cz = mri_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((mri_rrt_x-cx).^2 + (mri_rrt_y-cy).^2 + (mri_rrt_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
mri_rrt_dis = mri_rrt(any_inds==1); 

% STB mask velocity points (STB voxel averaged resolution)
any_inds = zeros(size(stb_rrt_x));
% Load STB velocity field
for q = 1:1:size(stb_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = stb_any_mask(q,1);
    cy = stb_any_mask(q,2);
    cz = stb_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((stb_rrt_x-cx).^2 + (stb_rrt_y-cy).^2 + (stb_rrt_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
stb_rrt_dis = stb_rrt(any_inds==1); 

% CFD mask velocity points (CFD voxel averaged resolution)
any_inds = zeros(size(cfd_rrt_x));
% Load STB velocity field
for q = 1:1:size(cfd_any_mask,1)
    % Get the current x, y, z points in the mask
    cx = cfd_any_mask(q,1);
    cy = cfd_any_mask(q,2);
    cz = cfd_any_mask(q,3);
    % Find the point in the current field closest to this
    pt_dist = sqrt((cfd_rrt_x-cx).^2 + (cfd_rrt_y-cy).^2 + (cfd_rrt_z-cz).^2);
    [~,close_ind] = min(pt_dist);
    any_inds(close_ind) = 1;
end
cfd_rrt_dis = cfd_rrt(any_inds==1);

% STB mask velocity points (STB full resolution)
% Using the velocity mask, identify all points whose closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(stb_maskX_va,stb_maskY_va,stb_maskZ_va);
stb_vamask = reshape(stb_velmask_va,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskX = reshape(Xva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskY = reshape(Yva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
stb_vamaskZ = reshape(Zva,[size(stb_velmask_va,1)*size(stb_velmask_va,2)*size(stb_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(stb_fres_rrt_x));
for q = 1:1:size(stb_fres_rrt_x)
    % Get current x,y,z point
    cx = stb_fres_rrt_x(q);
    cy = stb_fres_rrt_y(q);
    cz = stb_fres_rrt_z(q);
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
stb_rrt_fres_dis = stb_fres_rrt(fres_inany==1);


% CFD mask velocity points (CFD full resolution)
% Using the velocity mask, identify all points who's closest point on the
% voxel averaged field is in the aneurysm mask
% Start by identifying the voxel averaged mask points that are in the
% aneurysm
[Xva,Yva,Zva] = meshgrid(cfd_maskX_va,cfd_maskY_va,cfd_maskZ_va);
cfd_vamask = reshape(cfd_velmask_va,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskX = reshape(Xva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskY = reshape(Yva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
cfd_vamaskZ = reshape(Zva,[size(cfd_velmask_va,1)*size(cfd_velmask_va,2)*size(cfd_velmask_va,3),1]);
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
% Iterate through each point in the full resolution wss vector. Add full
% resolution point if its closest voxel averaged point is in the aneurysm
% mask
fres_inany = zeros(size(cfd_fres_rrt_x));
for q = 1:1:size(cfd_fres_rrt_x)
    % Get current x,y,z point
    cx = cfd_fres_rrt_x(q);
    cy = cfd_fres_rrt_y(q);
    cz = cfd_fres_rrt_z(q);
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
cfd_rrt_fres_dis = cfd_fres_rrt(fres_inany==1);

% Plot the rrt
fp = load('wss_cmap.mat');
wss_cmap = fp.cmap;
max_color = 80;
fsize = 20;
figure(479); hold off; scatter3(mri_rrt_x,mri_rrt_y,mri_rrt_z,10,mri_rrt,'filled')
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); %colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(480); hold off; scatter3(stb_rrt_x,stb_rrt_y,stb_rrt_z,10,stb_rrt,'filled')
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); %colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(481); hold off; scatter3(cfd_rrt_x,cfd_rrt_y,cfd_rrt_z,10,cfd_rrt,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); %colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(482); hold off; scatter3(stb_fres_rrt_x,stb_fres_rrt_y,stb_fres_rrt_z,10,stb_fres_rrt,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); %colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
figure(483); hold off; scatter3(cfd_fres_rrt_x,cfd_fres_rrt_y,cfd_fres_rrt_z,10,cfd_fres_rrt,'filled');
caxis([0 max_color]); set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); %colormap(gca,wss_cmap);
axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;
print(479,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_mri_rrt.pdf'])
print(480,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_rrt.pdf'])
print(481,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_rrt.pdf'])
print(482,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_fres_rrt.pdf'])
print(483,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_fres_rrt.pdf'])

% Compute time and space averaged TAWSS
mri_rrt_avg = mean(mri_rrt_dis);
stb_rrt_avg = mean(stb_rrt_dis);
cfd_rrt_avg = mean(cfd_rrt_dis);
stb_fres_rrt_avg = mean(stb_rrt_fres_dis);
cfd_fres_rrt_avg = mean(cfd_rrt_fres_dis);
mri_rrt_std = std(mri_rrt_dis);
stb_rrt_std = std(stb_rrt_dis);
cfd_rrt_std = std(cfd_rrt_dis);
stb_fres_rrt_std = std(stb_rrt_fres_dis);
cfd_fres_rrt_std = std(cfd_rrt_fres_dis);
rrt_avg = (1/10)*[mri_rrt_avg,stb_fres_rrt_avg,cfd_fres_rrt_avg,stb_rrt_avg,cfd_rrt_avg];
rrt_std = (1/10)*[mri_rrt_std,stb_fres_rrt_std,cfd_fres_rrt_std,stb_rrt_std,cfd_rrt_std];


%%% Get RRT statistics in aneurysm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_avg = [tawss_avg;osi_avg;rrt_avg];

%%% Create Violin Plots of data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAWSS violin plots
fsize = 22;
% Yviolin{1,1} = 10*tawss_all{1,1}; Yviolin{1,2} = 10*tawss_all{2,1};
% Yviolin{1,3} = 10*tawss_all{4,1}; Yviolin{1,4} = 10*tawss_all{5,1};
% Yviolin{1,5} = 10*tawss_all{3,1};
Yviolin{1,1} = 10*mri_tawss_dis; Yviolin{1,2} = 10*stb_tawss_dis;
Yviolin{1,3} = 10*stb_tawss_fres_dis; Yviolin{1,4} = 10*cfd_tawss_dis;
Yviolin{1,5} = 10*cfd_tawss_fres_dis;
%Yviolin = transpose(tawss_all);
[~,~,MX,MED,IQR,ci95,~] = violin(Yviolin,'figNumber',250,'edgecolor','none','facecolor',[0.5 0.5 0.5],'xlabel',{'4D Flow','STB-VA','STB','CFD-VA','CFD'},'plotlegend',0,'95CIdistribution','poisson');
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); ylabel('TAWSS (dynes/cm^{2})'); set(gca,'Box','off'); ylim([-2 40])
set(gca,'XTickLabel',{'4D Flow','STB-VA','STB','CFD-VA','CFD'}); pbaspect([1 0.6 1]);
ind1 = 1;
avg_all(ind1,:) = MX;
med_all(ind1,:) = MED;
IQR_all(1:2,:) = IQR;
ci95_all(1:2,:) = ci95;

% OSI violin plot
Yviolin{1,1} = mri_osi_dis; Yviolin{1,2} = stb_osi_dis;
Yviolin{1,3} = stb_osi_fres_dis; Yviolin{1,4} = cfd_osi_dis;
Yviolin{1,5} = cfd_osi_fres_dis;
[~,~,MX,MED,IQR,ci95,~] = violin(Yviolin,'figNumber',251,'edgecolor','none','facecolor',[0.5 0.5 0.5],'xlabel',{'4D Flow','STB-VA','STB','CFD-VA','CFD'},'plotlegend',0,'95CIdistribution','poisson');
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); ylabel('OSI'); set(gca,'Box','off');
if strcmp(cadat.FILE.anynum,'007'), ylim([-0.2 0.6]); 
else, ylim([-0.2 0.6]); end
ind1 = 2;
avg_all(ind1,:) = MX;
med_all(ind1,:) = MED;
IQR_all(3:4,:) = IQR;
ci95_all(3:4,:) = ci95;

set(gca,'XTickLabel',{'4D Flow','STB-VA','STB','CFD-VA','CFD'}); pbaspect([1 0.6 1]);
% RRT violin plot
% Yviolin{1,1} = (1/10)*rrt_all{1,1}; Yviolin{1,2} = (1/10)*rrt_all{2,1};
% Yviolin{1,3} = (1/10)*rrt_all{4,1}; Yviolin{1,4} = (1/10)*rrt_all{5,1};
% Yviolin{1,5} = (1/10)*rrt_all{3,1};
Yviolin{1,1} = (1/10)*mri_rrt_dis; Yviolin{1,2} = (1/10)*stb_rrt_dis;
Yviolin{1,3} = (1/10)*stb_rrt_fres_dis; Yviolin{1,4} = (1/10)*cfd_rrt_dis;
Yviolin{1,5} = (1/10)*cfd_rrt_fres_dis;
[~,~,MX,MED,IQR,ci95,~] = violin(Yviolin,'figNumber',252,'edgecolor','none','facecolor',[0.5 0.5 0.5],'xlabel',{'4D Flow','STB-VA','STB','CFD-VA','CFD'},'plotlegend',0,'95CIdistribution','poisson');
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); ylabel('RRT (dynes/cm^{2})^{-1}'); set(gca,'Box','off');
if strcmp(cadat.FILE.anynum,'111'), ylim([-3 10]); 
else, ylim([-3 10]); end
set(gca,'XTickLabel',{'4D Flow','STB-VA','STB','CFD-VA','CFD'}); pbaspect([1 0.6 1]);
ind1 = 3;
avg_all(ind1,:) = MX;
med_all(ind1,:) = MED;
IQR_all(5:6,:) = IQR;
ci95_all(5:6,:) = ci95;

print(250,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_TAWSS.pdf'])
print(251,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_OSI.pdf'])
print(252,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_RRT.pdf'])
savefig(250,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_TAWSS.fig'])
savefig(251,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_OSI.fig'])
savefig(252,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_RRT.fig'])

h = legend('PDF','Mean','Median','IQR','95% CI','Location','northoutside','Orientation','horizontal');
set(h,'Box','off');
print(252,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_ViolinPlot_legend.pdf'])

ww = 1;

% % % % Create bar plot of the time and space averaged TAWSS, OSI, RRT
% % % fsize = 24;
% % % val_avg = tawss_avg;
% % % bar_x = 0.4:0.4:2.6;
% % % bar_x = bar_x(1:length(val_avg));
% % % for q = 1:1:length(val_avg)
% % %     c_color = color_order(q,:);
% % %     figure(500); hold on; bar(bar_x(q),val_avg(q),'BarWidth',0.4,'FaceColor',c_color);
% % % end
% % % axis([0 2.4 0 20]); set(gcf,'Color',[1 1 1]); set(gca,'FontSize',fsize);
% % % set(gca,'XTick',bar_x)
% % % %set(gca,'XTickLabel',{'MRI','STB-Voxel','CFD-Voxel','STB-Full','CFD-Full'})
% % % set(gca,'XTickLabel',{}); set(gca,'Box','on'); grid on;
% % % ylabel('WSS (dyne/cm^{2})')
% % % pbaspect([0.7 1 1])
% % % 
% % % % OSI
% % % val_avg = osi_avg;
% % % bar_x = 0.4:0.4:3;
% % % bar_x = bar_x(1:length(val_avg));
% % % for q = 1:1:length(val_avg)
% % %     c_color = color_order(q,:);
% % %     figure(501); hold on; bar(bar_x(q),val_avg(q),'BarWidth',0.4,'FaceColor',c_color);
% % % end
% % % axis([0 2.4 0 0.15]); set(gcf,'Color',[1 1 1]); set(gca,'FontSize',fsize);
% % % set(gca,'XTick',bar_x)
% % % %set(gca,'XTickLabel',{'MRI','STB-Voxel','CFD-Voxel','STB-Full','CFD-Full'})
% % % set(gca,'XTickLabel',{}); set(gca,'Box','on'); grid on;
% % % ylabel('OSI')
% % % pbaspect([0.7 1 1])
% % % 
% % % % RRT
% % % val_avg = rrt_avg;
% % % bar_x = 0.4:0.4:3;
% % % bar_x = bar_x(1:length(val_avg));
% % % for q = 1:1:length(val_avg)
% % %     c_color = color_order(q,:);
% % %     figure(502); hold on; bar(bar_x(q),val_avg(q),'BarWidth',0.4,'FaceColor',c_color);
% % % end
% % % axis([0 2.4 0 2]); set(gcf,'Color',[1 1 1]); set(gca,'FontSize',fsize);
% % % set(gca,'XTick',bar_x)
% % % %set(gca,'XTickLabel',{'MRI','STB-Voxel','CFD-Voxel','STB-Full','CFD-Full'})
% % % set(gca,'XTickLabel',{}); set(gca,'Box','on'); grid on;
% % % ylabel('RRT (dyne/cm^{2})^{-1}')
% % % pbaspect([0.7 1 1])
% % % 
% % % print(500,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_barcomp_wss.pdf'])
% % % print(501,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_barcomp_osi.pdf'])
% % % print(502,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_barcomp_rrt.pdf'])
% % % 
% % % % Print the legend
% % % h = legend('4D Flow','STB-Full','CFD-Full','Orientation','horizontal','Location','northoutside');
% % % set(h,'Box','off')
% % % print(502,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_barcomp_legendFront.pdf'])
% % % h = legend('','','','STB-Voxel','CFD-Voxel','Orientation','horizontal','Location','northoutside');
% % % set(h,'Box','off')
% % % print(502,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_barcomp_legendBack.pdf'])




% Plot the histogram of TAWSS, OSI, and RRT
% Voxel average bins
mri_tawss_dis = 10*mri_tawss_dis;
stb_tawss_dis = 10*stb_tawss_dis;
cfd_tawss_dis = 10*cfd_tawss_dis;
stb_tawss_fres_dis = 10*stb_tawss_fres_dis;
cfd_tawss_fres_dis = 10*cfd_tawss_fres_dis;
num_bins = round(length(mri_tawss_dis)/10);%8;
min_bin_va = min([min(cfd_tawss_dis),min(stb_tawss_dis),min(mri_tawss_dis),min(stb_tawss_fres_dis),min(cfd_tawss_fres_dis)]);
max_bin_va = max([max(cfd_tawss_dis),max(stb_tawss_dis),max(mri_tawss_dis),max(stb_tawss_fres_dis),max(cfd_tawss_fres_dis)]);
last_bin_edge = 30; %ceil(max_bin_va);
bin_edges_va = linspace(0,last_bin_edge,num_bins+1);
%bin_edges_va = [bin_edges_va,max_bin_va];
plt_bins_va = (bin_edges_va(2:end)-bin_edges_va(1:end-1))/2 + bin_edges_va(1:end-1);
% Full resolution bins
num_bins = round(length(cfd_tawss_fres_dis)/50);%25;
min_bin = min([min(cfd_tawss_dis),min(stb_tawss_dis),min(stb_tawss_fres_dis),min(cfd_tawss_fres_dis),min(mri_tawss_dis)]);
max_bin = max([max(cfd_tawss_dis),max(stb_tawss_dis),max(stb_tawss_fres_dis),max(cfd_tawss_fres_dis),max(mri_tawss_dis)]);
last_bin_edge = 30; %ceil(max_bin);
bin_edges = linspace(0,last_bin_edge,num_bins+1);
%bin_edges = [bin_edges,max_bin];
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
% Compute histogram
figure(1);
stb_pdf = histogram(stb_tawss_dis,bin_edges_va,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_tawss_dis,bin_edges_va,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_tawss_fres_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_pdf = histogram(cfd_tawss_dis,bin_edges_va,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_tawss_dis,bin_edges_va,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
cfd_fres_pdf = histogram(cfd_tawss_fres_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
mri_pdf = histogram(mri_tawss_dis,bin_edges_va,'Normalization','pdf');
mri_pdf = mri_pdf.Values;



% Plot the PDFS
lsize = 3;
fsize = 24;
max_val = max([mri_pdf,cfd_pdf,stb_pdf,cfd_fres_pdf,stb_fres_pdf]);
figure(611); hold off; plot(plt_bins_va,mri_pdf/max_val,'Color',cmri,'LineWidth',lsize);
hold on; plot(plt_bins_va,stb_pdf/max_val,'--','Color',cstb,'LineWidth',lsize);
plot(plt_bins_va,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',lsize);
plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',lsize);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',lsize);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); set(gca,'Box','off')
xlabel('TAWSS (dynes/cm^{2})')
ylabel('PDF')
pbaspect([1 0.6 1])
print(611,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_tawss_pdf.pdf'])


% Voxel average bins
num_bins = round(length(mri_osi_dis)/10);%8;
min_bin_va = min([min(cfd_osi_dis),min(stb_osi_dis),min(mri_osi_dis),min(stb_osi_fres_dis),min(cfd_osi_fres_dis)]);
max_bin_va = max([max(cfd_osi_dis),max(stb_osi_dis),max(mri_osi_dis),max(stb_osi_fres_dis),max(cfd_osi_fres_dis)]);
last_bin_edge = 0.5; %ceil(max_bin_va);
bin_edges_va = linspace(0,last_bin_edge,num_bins+1);
%bin_edges_va = [bin_edges_va,max_bin_va];
plt_bins_va = (bin_edges_va(2:end)-bin_edges_va(1:end-1))/2 + bin_edges_va(1:end-1);
% Full resolution bins
num_bins = round(length(cfd_osi_fres_dis)/50);%25;
min_bin = min([min(cfd_osi_dis),min(stb_osi_dis),min(stb_osi_fres_dis),min(cfd_osi_fres_dis),min(mri_osi_dis)]);
max_bin = max([max(cfd_osi_dis),max(stb_osi_dis),max(stb_osi_fres_dis),max(cfd_osi_fres_dis),max(mri_osi_dis)]);
last_bin_edge = 0.5; %ceil(max_bin);
bin_edges = linspace(0,last_bin_edge,num_bins+1);
%bin_edges = [bin_edges,max_bin];
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
% Compute histogram
figure(1);
stb_pdf = histogram(stb_osi_dis,bin_edges_va,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_cdf = histogram(stb_osi_dis,bin_edges_va,'Normalization','cdf');
stb_cdf = stb_cdf.Values;
stb_fres_pdf = histogram(stb_osi_fres_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_pdf = histogram(cfd_osi_dis,bin_edges_va,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_cdf = histogram(cfd_osi_dis,bin_edges_va,'Normalization','cdf');
cfd_cdf = cfd_cdf.Values;
cfd_fres_pdf = histogram(cfd_osi_fres_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
mri_pdf = histogram(mri_osi_dis,bin_edges_va,'Normalization','pdf');
mri_pdf = mri_pdf.Values;

% Plot the PDFS
lsize = 3;
fsize = 24;
max_val = max([mri_pdf,cfd_pdf,stb_pdf,cfd_fres_pdf,stb_fres_pdf]);
figure(612); hold off; plot(plt_bins_va,mri_pdf/max_val,'Color',cmri,'LineWidth',lsize);
hold on; plot(plt_bins_va,stb_pdf/max_val,'--','Color',cstb,'LineWidth',lsize);
plot(plt_bins_va,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',lsize);
plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',lsize);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',lsize);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); set(gca,'Box','off')
xlabel('OSI')
ylabel('PDF')
pbaspect([1 0.6 1])
print(612,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_osi_pdf.pdf'])


% Voxel average bins
mri_rrt_dis = (1/10)*mri_rrt_dis;
stb_rrt_dis = (1/10)*stb_rrt_dis;
cfd_rrt_dis = (1/10)*cfd_rrt_dis;
stb_rrt_fres_dis = (1/10)*stb_rrt_fres_dis;
cfd_rrt_fres_dis = (1/10)*cfd_rrt_fres_dis;
num_bins = round(length(mri_rrt_dis)/10);%8;
min_bin_va = min([min(cfd_rrt_dis),min(stb_rrt_dis),min(mri_rrt_dis),min(stb_rrt_fres_dis),min(cfd_rrt_fres_dis)]);
max_bin_va = max([max(cfd_rrt_dis),max(stb_rrt_dis),max(mri_rrt_dis),max(stb_rrt_fres_dis),max(cfd_rrt_fres_dis)]);
last_bin_edge = 5; %ceil(max_bin_va);
bin_edges_va = linspace(0,last_bin_edge,num_bins+1);
%bin_edges_va = [bin_edges_va,max_bin_va];
plt_bins_va = (bin_edges_va(2:end)-bin_edges_va(1:end-1))/2 + bin_edges_va(1:end-1);
% Full resolution bins
num_bins = round(length(cfd_rrt_fres_dis)/50);%25;
min_bin = min([min(cfd_rrt_dis),min(stb_rrt_dis),min(stb_rrt_fres_dis),min(cfd_rrt_fres_dis),min(mri_rrt_dis)]);
max_bin = max([max(cfd_rrt_dis),max(stb_rrt_dis),max(stb_rrt_fres_dis),max(cfd_rrt_fres_dis),max(mri_rrt_dis)]);
last_bin_edge = 5; %ceil(max_bin);
bin_edges = linspace(0,last_bin_edge,num_bins+1);
%bin_edges = [bin_edges,max_bin];
plt_bins = (bin_edges(2:end)-bin_edges(1:end-1))/2 + bin_edges(1:end-1);
% Compute histogram
figure(1);
stb_pdf = histogram(stb_rrt_dis,bin_edges_va,'Normalization','pdf');
stb_pdf = stb_pdf.Values;
stb_fres_pdf = histogram(stb_rrt_fres_dis,bin_edges,'Normalization','pdf');
stb_fres_pdf = stb_fres_pdf.Values;
cfd_pdf = histogram(cfd_rrt_dis,bin_edges_va,'Normalization','pdf');
cfd_pdf = cfd_pdf.Values;
cfd_fres_pdf = histogram(cfd_rrt_fres_dis,bin_edges,'Normalization','pdf');
cfd_fres_pdf = cfd_fres_pdf.Values;
mri_pdf = histogram(mri_rrt_dis,bin_edges_va,'Normalization','pdf');
mri_pdf = mri_pdf.Values;

% Plot the PDFS
lsize = 3;
fsize = 24;
max_val = max([mri_pdf,cfd_pdf,stb_pdf,cfd_fres_pdf,stb_fres_pdf]);
figure(613); hold off; plot(plt_bins_va,mri_pdf/max_val,'Color',cmri,'LineWidth',lsize);
hold on; plot(plt_bins_va,stb_pdf/max_val,'--','Color',cstb,'LineWidth',lsize);
plot(plt_bins_va,cfd_pdf/max_val,'--','Color',ccfd,'LineWidth',lsize);
plot(plt_bins,stb_fres_pdf/max_val,'Color',cstb_light,'LineWidth',lsize);
plot(plt_bins,cfd_fres_pdf/max_val,'Color',ccfd_light,'LineWidth',lsize);
set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); set(gca,'Box','off')
xlabel('RRT (dynes/cm^{2})^{-1}')
ylabel('PDF')
pbaspect([1 0.6 1])
print(613,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_any_rrt_pdf.pdf'])



%%% BLAND-ALTMAN OF TAWSS, OSI, RRT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This analysis is done only comparing STB/CFD voxel vs full resolution
%%% PLOT BLAND-ALTMAN COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plotting is done at peak systole, peak diastole, and mid diastole in 
% order to compare all modalities
%mean_color = [135,206,250]/255;
mean_color = [65,105,225]/255;
mean_color_2 = [65,161,225]/255;
% Iterate through the three TAWSS, OSI, RRT
% % % for q = 1:1:3
% % %     % Set the parameters based on the q
% % %     if q == 1
% % %         u_all = tawss_all;
% % %         x_all = tawss_all_X;
% % %         y_all = tawss_all_Y;
% % %         z_all = tawss_all_Z;
% % %     elseif q == 2
% % %         u_all = osi_all;
% % %         x_all = osi_all_X;
% % %         y_all = osi_all_Y;
% % %         z_all = osi_all_Z;
% % %     else
% % %         u_all = rrt_all;
% % %         x_all = rrt_all_X;
% % %         y_all = rrt_all_Y;
% % %         z_all = rrt_all_Z;
% % %     end
% % %     % Iterate through all file types
% % %     iter = 1; % Set the iteration counter
% % %     if strcmp(cadat.FILE.anynum,'111')
% % %         % Location comparison limits
% % %         minX = -10.15; maxX = 7.55;
% % %         minY = -5.15; maxY = 13.75;
% % %         minZ = -7; maxZ = 7;
% % %     else
% % %         minX = -14; maxX = 8;
% % %         minY = -6; maxY = 24;
% % %         minZ = -10; maxZ = 15;
% % %     end
% % %     % Subplot location settings
% % %     dc = 0.19; c2 = 0.5025; %0.13+dc; c3 = 0.13+2*dc;
% % %     dr = 0.25; r2 = 0.5025; %0.7093-dr; r3 = 0.7093-2*dr;
% % %     %umin = -0.25; umax = 0.25; % Max and min for u-vel x-axis settings
% % %     %wmin = -0.15; wmax = 0.15; % Max and min for u-vel x-axis settings
% % %     %ymin = -0.25; ymax = 0.25; % Max and min for y-axis settings
% % %     if strcmp(cadat.FILE.anynum,'007')
% % %         if q == 1
% % %             scale_fact = 10;
% % %             ymin = -20; ymax = 5; % Max and min for y-axis settings
% % %             dx = 3;
% % %             x_axis_label = 'mean (dynes/cm^{2})';
% % %             y_axis_label = 'difference (dynes/cm^{2})';
% % %         elseif q == 2
% % %             scale_fact = 1;
% % %             ymin = -0.4; ymax = 0.4; % Max and min for y-axis settings
% % %             dx = 0.1;
% % %             x_axis_label = 'mean';
% % %             y_axis_label = 'difference';
% % %         else
% % %             scale_fact = 0.1;
% % %             ymin = -40; ymax = 40; % Max and min for y-axis settings
% % %             dx = 3;
% % %             x_axis_label = 'mean (dynes/cm^{2})^{-1}';
% % %             y_axis_label = 'difference (dynes/cm^{2})^{-1}';
% % %         end
% % %         xmax = [12,0.5,12];
% % %     else
% % %         if q == 1
% % %             scale_fact = 10;
% % %             ymin = -30; ymax = 4; % Max and min for y-axis settings
% % %             dx = 6;
% % %             x_axis_label = 'mean (dynes/cm^{2})';
% % %             y_axis_label = 'difference (dynes/cm^{2})';
% % %         elseif q == 2
% % %             scale_fact = 1;
% % %             ymin = -0.4; ymax = 0.4; % Max and min for y-axis settings
% % %             dx = 0.1;
% % %             x_axis_label = 'mean';
% % %             y_axis_label = 'difference';
% % %         else
% % %             scale_fact = 0.1;
% % %             ymin = -20; ymax = 20; % Max and min for y-axis settings
% % %             dx = 2;
% % %             x_axis_label = 'mean (dynes/cm^{2})^{-1}';
% % %             y_axis_label = 'difference (dynes/cm^{2})^{-1}';
% % %         end
% % %         xmax = [24,0.5,8];
% % %     end
% % %     
% % %     %umin = -0.4; umax = 0.4; % Max and min for u-vel x-axis settings
% % %     for ftype1 = 2:1:3
% % %         % Set the data for the first file type
% % %         x1 = x_all(ftype1,1); x1 = x1{1,1};
% % %         y1 = y_all(ftype1,1); y1 = y1{1,1};
% % %         z1 = z_all(ftype1,1); z1 = z1{1,1};
% % %         u1 = u_all(ftype1,1); u1 = u1{1,1};
% % %         
% % %         % Set the aneurysm mask points for the first time point
% % %         if ftype1 == 1
% % %             % Set the aneurysm mask for the first dataset
% % %             any_mask_1 = mri_any_mask;
% % %             % Set the velocity mask information for the current dataset
% % %             xmask_1 = mri_maskX;
% % %             ymask_1 = mri_maskY;
% % %             zmask_1 = mri_maskZ;
% % %         elseif ftype1 == 2
% % %             any_mask_1 = stb_any_mask;
% % %             % Set the velocity mask information for the current dataset
% % %             xmask_1 = stb_maskX;
% % %             ymask_1 = stb_maskY;
% % %             zmask_1 = stb_maskZ;
% % %         else
% % %             any_mask_1 = cfd_any_mask;
% % %             % Set the velocity mask information for the current dataset
% % %             xmask_1 = cfd_maskX;
% % %             ymask_1 = cfd_maskY;
% % %             zmask_1 = cfd_maskZ;
% % %         end
% % %         
% % %         % Iterate through all other modalities to compare
% % %         for ftype2 = ftype1+2:1:ftype1+2
% % %             % Set the minimum threshold
% % %             if ftype1 == 2 || ftype2 == 2 || ftype1 == 3 || ftype2 == 3
% % %                 grid_half = 0.67; %0.15;
% % %                 dist_thresh = sqrt(3*grid_half^2);
% % %             else
% % %                 grid_half = 0.67;
% % %                 dist_thresh = sqrt(3*grid_half^2);
% % %             end
% % %             
% % %             % Set the aneurysm mask points for the first time point
% % %             if ftype2 == 4
% % %                 % Set the aneurysm mask for the first dataset
% % %                 any_mask_2 = stb_any_mask;
% % %                 % Set the velocity mask information for the current dataset
% % %                 xmask_2 = stb_maskX;
% % %                 ymask_2 = stb_maskY;
% % %                 zmask_2 = stb_maskZ;
% % %             else
% % %                 any_mask_2 = cfd_any_mask;
% % %                 % Set the velocity mask information for the current dataset
% % %                 xmask_2 = cfd_maskX;
% % %                 ymask_2 = cfd_maskY;
% % %                 zmask_2 = cfd_maskZ;
% % %             end        
% % %             
% % %             % Set the data for the first file type
% % %             x2 = x_all(ftype2,1); x2 = x2{1,1};
% % %             y2 = y_all(ftype2,1); y2 = y2{1,1};
% % %             z2 = z_all(ftype2,1); z2 = z2{1,1};
% % %             u2 = u_all(ftype2,1); u2 = u2{1,1};
% % %             
% % %             % Set the primary array as the shorter one, secondary as longer
% % %             % Set the aneurysm mask to be used based on the primary dataset
% % %             if length(x1)>length(x2)
% % %                 xp = x2; yp = y2; zp = z2; up = u2; 
% % %                 xs = x1; ys = y1; zs = z1; us = u1;
% % %                 any_mask = any_mask_2;
% % %                 xmask = xmask_2; ymask = ymask_2; zmask = zmask_2;
% % %                 mult_factor = -1;
% % %             else
% % %                 xs = x2; ys = y2; zs = z2; us = u2;
% % %                 xp = x1; yp = y1; zp = z1; up = u1;
% % %                 any_mask = any_mask_1;
% % %                 xmask = xmask_1; ymask = ymask_1; zmask = zmask_1;
% % %                 mult_factor = 1;
% % %             end
% % %             
% % %             %%% Iterate through all data points to create the B-A arrays
% % %             ct = 1;
% % %             % Initialize outputs
% % %             dist_disparity = [];
% % %             ba_vec1 = [];
% % %             ba_vec2 = [];
% % %             for pn = 1:1:length(xp)
% % %                 % Get the current point from the x primary
% % %                 cx = xp(pn);
% % %                 cy = yp(pn);
% % %                 cz = zp(pn);
% % %                 
% % %                 % Compute the distances from the current point
% % %                 pt_dists = sqrt((xs-cx).^2 + (ys-cy).^2 + (zs-cz).^2);
% % %                 [min_dist,sn] = min(pt_dists); % Get the minimum distance
% % %                 
% % %                 % Determine if the current value is in range
% % %                 in_range_x = (cx >= minX) & (cx <= maxX);
% % %                 in_range_y = (cy >= minY) & (cy <= maxY);
% % %                 in_range_z = (cz >= minZ) & (cz <= maxZ);
% % %                 all_in_range = in_range_x & in_range_y & in_range_z;
% % %                 
% % %                 % Determine if the current point is in the aneurysm sac
% % %                 % Find the closest point in the mask to the current point
% % %                 [~,min_ptX] = min(abs(xmask-cx)); % Find the closest mask point
% % %                 [~,min_ptY] = min(abs(ymask-cy)); % Find the closest mask point
% % %                 [~,min_ptZ] = min(abs(zmask-cz)); % Find the closest mask point
% % %                 xm_pt = xmask(min_ptX); 
% % %                 ym_pt = ymask(min_ptY);
% % %                 zm_pt = zmask(min_ptZ);
% % %                 % Determine if the current mask point is in the aneurysm
% % %                 % mask array
% % %                 m_pt = repmat([xm_pt,ym_pt,zm_pt],[size(any_mask,1),1]);
% % %                 m_pt_dist = sum(abs(m_pt-any_mask),2);
% % %                 min_pt_dist = min(m_pt_dist);
% % %                 in_any_mask = min_pt_dist < grid_half; % Determine if the point is in the aneurysm mask
% % %                 
% % %                 % If the minimum distance is below the threshold and the 
% % %                 % point is within the defined ROI, a match is accepted, add
% % %                 % the point to the output arrays
% % %                 if (min_dist < dist_thresh) && all_in_range && in_any_mask
% % %                     dist_disparity(ct,1) = min_dist;
% % %                     ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
% % %                     ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
% % %                     ct = ct + 1;
% % %                 end
% % %             end
% % %             
% % %             % Set index skip for plotting
% % %             iskip = 1; %if iter < 3, iskip = 1;
% % %             %else iskip = 10; end
% % %             
% % %             % Compute the B-A variables and  plot them
% % %             pbrat = 1; % Size ratio of plots
% % %             msize = 4; % Set the marker size
% % %             line_size = 2; % Set the line width size
% % %             font_size = 16; % Set the font size
% % % 
% % %             %%% plot B-A analysis
% % %             ba_vec1(:,4) = scale_fact*ba_vec1(:,4);
% % %             ba_vec2(:,4) = scale_fact*ba_vec2(:,4);
% % %             ba_mean = mean([ba_vec1(:,4),ba_vec2(:,4)],2);
% % %             ba_diff = mult_factor*(ba_vec1(:,4) - ba_vec2(:,4));
% % %             %ba_diff = abs(ba_diff);
% % %             %ba_diff = ba_diff_act./ba_mean;
% % %             meanDiff = mean(ba_diff);
% % %             sdDiff = std(ba_diff);
% % %             CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
% % %             % Compute linear fit of means and differences
% % %             linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
% % %             % Plot the u-velocity B-A analysis
% % %             figure(600+q); h1 = subplot(1,2,iter);
% % %             hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
% % %             plot([-20 0 60], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
% % %             plot([-20 0 60], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
% % %             %plot([-20 20], [-20 20].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
% % %             plot([-20 0 60], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
% % %             axis([0 xmax(q) ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
% % %             set(gca,'XTick',0:dx:xmax(q)); grid('on');
% % %             if iter == 2
% % %                 set(gca,'YTickLabel',{});
% % %             else
% % %                 ylabel(y_axis_label)
% % %             end
% % %             % Set the axis information
% % %             xlabel(x_axis_label)
% % %             
% % %                 
% % %             p = get(h1,'Position');
% % %             % Set the subplot location (row)
% % %             if iter > 1, p(1) = c2; end
% % %             %elseif 3*(iter-1)+1 > 3, p(2) = r2;
% % %             %end
% % %             set(h1,'Position',p);
% % %             % Set title
% % %             if q == 1, title('TAWSS'); 
% % %             elseif q == 2, title('OSI');
% % %             else, title('RRT'); 
% % %             end
% % %             % Save info in array
% % %             ba_results{q,1}(:,iter) = [meanDiff;CR(1);CR(2)];
% % %             
% % %             
% % %             % Save the disparity error in an array
% % %             dist_disparity_error{iter,1} = dist_disparity;
% % %             dist_disparity_avgerror(iter,1) = mean(dist_disparity);        
% % %             
% % %             % Increment the iteration counter
% % %             iter = iter + 1;    
% % %         end
% % %     end
% % %     % Get the cycle point print info
% % %     if q == 1, cycle_pt = 'TAWSS_';
% % %     elseif q == 2, cycle_pt = 'OSI_';
% % %     else, cycle_pt = 'RRT_'; 
% % %     end   
% % %     % Print the Bland-Altman plots
% % %     print(600+q,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-anymask-comp.pdf'])
% % % end


%%% BLAND-ATLMAN USING THE ENTIRE DOMAIN, NOT JUST THE ANEURYSM SAC %%%%%%%
for q = 1:1:3
    % Set the parameters based on the q
    if q == 1
        u_all = tawss_all;
        x_all = tawss_all_X;
        y_all = tawss_all_Y;
        z_all = tawss_all_Z;
    elseif q == 2
        u_all = osi_all;
        x_all = osi_all_X;
        y_all = osi_all_Y;
        z_all = osi_all_Z;
    else
        u_all = rrt_all;
        x_all = rrt_all_X;
        y_all = rrt_all_Y;
        z_all = rrt_all_Z;
    end
    
    iter = 1; % Set the iteration counter
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
    % Subplot location settings
    dc = 0.19; c2 = 0.5025; %0.13+dc; c3 = 0.13+2*dc;
    dr = 0.25; r2 = 0.5025; %0.7093-dr; r3 = 0.7093-2*dr;
    %umin = -0.25; umax = 0.25; % Max and min for u-vel x-axis settings
    %wmin = -0.15; wmax = 0.15; % Max and min for u-vel x-axis settings
    %ymin = -0.25; ymax = 0.25; % Max and min for y-axis settings
    if q == 1
        scale_fact = 10;
        ymin = -30; ymax = 10; % Max and min for y-axis settings
        dx = 6;
        x_axis_label = 'mean (dynes/cm^{2})';
        y_axis_label = 'difference (dynes/cm^{2})';
    elseif q == 2
        scale_fact = 1;
        ymin = -0.4; ymax = 0.4; % Max and min for y-axis settings
        dx = 0.2;
        x_axis_label = 'mean';
        y_axis_label = 'difference';
    else
        scale_fact = 0.1;
        ymin = -12; ymax = 30; % Max and min for y-axis settings
        dx = 5;
        x_axis_label = 'mean (dynes/cm^{2})^{-1}';
        y_axis_label = 'difference (dynes/cm^{2})^{-1}';
    end
    xmax = [24,0.5,15];
    
    for ftype1 = 2:1:3
        % Set the data for the first file type
        x1 = x_all(ftype1,1); x1 = x1{1,1};
        y1 = y_all(ftype1,1); y1 = y1{1,1};
        z1 = z_all(ftype1,1); z1 = z1{1,1};
        u1 = u_all(ftype1,1); u1 = u1{1,1};
        % Iterate through all other modalities to compare
        for ftype2 = ftype1+2:1:ftype1+2
            % Set the minimum threshold
            if ftype1 == 2 || ftype2 == 2 || ftype1 == 3 || ftype2 == 3
                grid_half = 0.67; %%0.15;
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

            % Set the primary array as the shorter one, secondary as longer
            if length(x1)>length(x2)
                xp = x2; yp = y2; zp = z2; up = u2;
                xs = x1; ys = y1; zs = z1; us = u1;
                mult_factor = -1;
            else
                xs = x2; ys = y2; zs = z2; us = u2; 
                xp = x1; yp = y1; zp = z1; up = u1;
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
                % If the minimum distance is below the threshold and the 
                % point is within the defined ROI, a match is accepted, add
                % the point to the output arrays
                if (min_dist < dist_thresh) && all_in_range
                    if q == 3
                        if up(pn) < 800 && us(sn) < 800
                            dist_disparity(ct,1) = min_dist;
                            ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
                            ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
                            ct = ct + 1;
                        end
                    else
                        dist_disparity(ct,1) = min_dist;
                        ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
                        ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
                        ct = ct + 1;
                    end
                end
            end

            % Set index skip for plotting
            iskip = 1; %if iter < 3, iskip = 1;
            %else iskip = 10; end

            % Compute the B-A variables and  plot them
            pbrat = 1; % Size ratio of plots
            msize = 2; % Set the marker size
            line_size = 2; % Set the line width size
            font_size = 16; % Set the font size

            %%% plot B-A analysis
            ba_vec1(:,4) = scale_fact*ba_vec1(:,4);
            ba_vec2(:,4) = scale_fact*ba_vec2(:,4);
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
            figure(200+q); h1 = subplot(1,2,iter);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-20 0 100], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-20 0 100], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-20 20], [-20 20].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-20 0 100], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
            %text(0.7*xmax(q),1.15*CR(1),['+CI = ',num2str(CR(1),'%.2f')]);
            %text(0.7*xmax(q),1.15*CR(2),['-CI = ',num2str(CR(2),'%.2f')]);
            %text(0.7*xmax(q),max((1.2*meanDiff),meanDiff+0.03),['MD = ',num2str(meanDiff,'%.2f')]);
            axis([0 xmax(q) ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'XTick',0:dx:xmax(q)); grid('on');
            if iter == 2
                set(gca,'YTickLabel',{});
            else
                ylabel(y_axis_label)
            end
            % Set the axis information
            xlabel(x_axis_label)


            p = get(h1,'Position');
            % Set the subplot location (row)
            if iter > 1, p(1) = c2; end
            %elseif 3*(iter-1)+1 > 3, p(2) = r2;
            %end
            set(h1,'Position',p);
            % Set title
            if q == 1, title('TAWSS'); 
            elseif q == 2, title('OSI');
            else, title('RRT'); 
            end
            % Save info in array
            ba_results{q,1}(:,iter) = [meanDiff;CR(1);CR(2)];


            % Save the disparity error in an array
            dist_disparity_error{iter,1} = dist_disparity;
            dist_disparity_avgerror(iter,1) = mean(dist_disparity);        

            % Increment the iteration counter
            iter = iter + 1;    
        end
    end
    % Get the cycle point print info
    if q == 1, cycle_pt = 'TAWSS_';
    elseif q == 2, cycle_pt = 'OSI_';
    else, cycle_pt = 'RRT_'; 
    end   
    % Print the Bland-Altman plots
    print(200+q,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-comp.pdf'])
end


%%% BLAND-ATLMAN USING THE ENTIRE DOMAIN, NOT JUST THE ANEURYSM SAC %%%%%%%
% % % for q = 1:1:3
% % %     % Set the parameters based on the q
% % %     if q == 1
% % %         u_all = tawss_all;
% % %         x_all = tawss_all_X;
% % %         y_all = tawss_all_Y;
% % %         z_all = tawss_all_Z;
% % %     elseif q == 2
% % %         u_all = osi_all;
% % %         x_all = osi_all_X;
% % %         y_all = osi_all_Y;
% % %         z_all = osi_all_Z;
% % %     else
% % %         u_all = rrt_all;
% % %         x_all = rrt_all_X;
% % %         y_all = rrt_all_Y;
% % %         z_all = rrt_all_Z;
% % %     end
% % %     
% % %     iter = 1; % Set the iteration counter
% % %     % Subplot location settings
% % %     dc = 0.1847; c2 = 0.13+dc; c3 = 0.13+2*dc; c4 = 0.13+3*dc;
% % %     dr = 0.25; r2 = 0.5025; %0.7093-dr; r3 = 0.7093-2*dr;
% % %     %umin = -0.25; umax = 0.25; % Max and min for u-vel x-axis settings
% % %     %wmin = -0.15; wmax = 0.15; % Max and min for u-vel x-axis settings
% % %     %ymin = -0.25; ymax = 0.25; % Max and min for y-axis settings
% % %     if q == 1
% % %         scale_fact = 10;
% % %         ymin = -32; ymax = 20; % Max and min for y-axis settings
% % %         dx = 6;
% % %         x_axis_label = 'mean (dynes/cm^{2})';
% % %         y_axis_label = 'difference (dynes/cm^{2})';
% % %     elseif q == 2
% % %         scale_fact = 1;
% % %         ymin = -0.4; ymax = 0.405; % Max and min for y-axis settings
% % %         dx = 0.2;
% % %         x_axis_label = 'mean';
% % %         y_axis_label = 'difference';
% % %     else
% % %         scale_fact = 0.1;
% % %         ymin = -30; ymax = 30; % Max and min for y-axis settings
% % %         dx = 2;
% % %         x_axis_label = 'mean (dynes/cm^{2})^{-1}';
% % %         y_axis_label = 'difference (dynes/cm^{2})^{-1}';
% % %     end
% % %     xmax = [24,0.4,10];
% % %     
% % %     for ftype1 = 1
% % %         % Set the data for the first file type
% % %         x1 = x_all(ftype1,1); x1 = x1{1,1};
% % %         y1 = y_all(ftype1,1); y1 = y1{1,1};
% % %         z1 = z_all(ftype1,1); z1 = z1{1,1};
% % %         u1 = u_all(ftype1,1); u1 = u1{1,1};
% % %         % Iterate through all other modalities to compare
% % %         for ftype2 = [2,4,3,5]
% % %             % Set the minimum threshold
% % %             if ftype1 == 2 || ftype2 == 2 || ftype1 == 3 || ftype2 == 3
% % %                 grid_half = 0.67; %%0.15;
% % %                 dist_thresh = sqrt(3*grid_half^2);
% % %             else
% % %                 grid_half = 0.67;
% % %                 dist_thresh = sqrt(3*grid_half^2);
% % %             end
% % % 
% % %             % Set the data for the first file type
% % %             x2 = x_all(ftype2,1); x2 = x2{1,1};
% % %             y2 = y_all(ftype2,1); y2 = y2{1,1};
% % %             z2 = z_all(ftype2,1); z2 = z2{1,1};
% % %             u2 = u_all(ftype2,1); u2 = u2{1,1};
% % % 
% % %             % Set the primary array as the shorter one, secondary as longer
% % %             if length(x1)>length(x2)
% % %                 xp = x2; yp = y2; zp = z2; up = u2;
% % %                 xs = x1; ys = y1; zs = z1; us = u1;
% % %             else
% % %                 xs = x2; ys = y2; zs = z2; us = u2; 
% % %                 xp = x1; yp = y1; zp = z1; up = u1;
% % %             end
% % % 
% % %             %%% Iterate through all data points to create the B-A arrays
% % %             ct = 1;
% % %             % Initialize outputs
% % %             dist_disparity = [];
% % %             ba_vec1 = [];
% % %             ba_vec2 = [];
% % %             for pn = 1:1:length(xp)
% % %                 % Get the current point from the x primary
% % %                 cx = xp(pn);
% % %                 cy = yp(pn);
% % %                 cz = zp(pn);
% % % 
% % %                 % Compute the distances from the current point
% % %                 pt_dists = sqrt((xs-cx).^2 + (ys-cy).^2 + (zs-cz).^2);
% % %                 [min_dist,sn] = min(pt_dists); % Get the minimum distance
% % % 
% % %                 % Determine if the current value is in range
% % %                 in_range_x = (cx >= minX) & (cx <= maxX);
% % %                 in_range_y = (cy >= minY) & (cy <= maxY);
% % %                 in_range_z = (cz >= minZ) & (cz <= maxZ);
% % %                 all_in_range = in_range_x & in_range_y & in_range_z;
% % %                 % If the minimum distance is below the threshold and the 
% % %                 % point is within the defined ROI, a match is accepted, add
% % %                 % the point to the output arrays
% % %                 if (min_dist < dist_thresh) && all_in_range
% % %                     if q == 3
% % %                         if up(pn) < 800 && us(sn) < 800
% % %                             dist_disparity(ct,1) = min_dist;
% % %                             ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
% % %                             ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
% % %                             ct = ct + 1;
% % %                         end
% % %                     else
% % %                         dist_disparity(ct,1) = min_dist;
% % %                         ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
% % %                         ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
% % %                         ct = ct + 1;
% % %                     end
% % %                 end
% % %             end
% % % 
% % %             % Set index skip for plotting
% % %             iskip = 1; %if iter < 3, iskip = 1;
% % %             %else iskip = 10; end
% % % 
% % %             % Compute the B-A variables and  plot them
% % %             pbrat = 1; % Size ratio of plots
% % %             msize = 1.5; % Set the marker size
% % %             line_size = 2; % Set the line width size
% % %             font_size = 10; % Set the font size
% % % 
% % %             %%% plot B-A analysis
% % %             ba_vec1(:,4) = scale_fact*ba_vec1(:,4);
% % %             ba_vec2(:,4) = scale_fact*ba_vec2(:,4);
% % %             ba_mean = mean([ba_vec1(:,4),ba_vec2(:,4)],2);
% % %             ba_diff = mult_factor*(ba_vec1(:,4) - ba_vec2(:,4));
% % %             %ba_diff = abs(ba_diff);
% % %             %ba_diff = ba_diff_act./ba_mean;
% % %             meanDiff = mean(ba_diff);
% % %             sdDiff = std(ba_diff);
% % %             CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff];
% % %             % Compute linear fit of means and differences
% % %             linFit = polyfit(ba_mean,ba_diff,1); %%%work out the linear fit coefficients
% % %             % Plot the u-velocity B-A analysis
% % %             figure(300+q); h1 = subplot(1,4,iter);
% % %             hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
% % %             plot([-20 0 100], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
% % %             plot([-20 0 100], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
% % %             %plot([-20 20], [-20 20].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
% % %             plot([-20 0 100], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
% % %             %text(0.7*xmax(q),1.15*CR(1),['+CI = ',num2str(CR(1),'%.2f')]);
% % %             %text(0.7*xmax(q),1.15*CR(2),['-CI = ',num2str(CR(2),'%.2f')]);
% % %             %text(0.7*xmax(q),max((1.2*meanDiff),meanDiff+0.03),['MD = ',num2str(meanDiff,'%.2f')]);
% % %             axis([0 xmax(q) ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
% % %             set(gca,'XTick',0:dx:xmax(q)); grid('on');
% % %             if iter > 1
% % %                 set(gca,'YTickLabel',{});
% % %             else
% % %                 ylabel(y_axis_label)
% % %             end
% % %             % Set the axis information
% % %             xlabel(x_axis_label)
% % % 
% % % 
% % %             p = get(h1,'Position');
% % %             % Set the subplot location (row)
% % %             if iter == 2, p(1) = c2;
% % %             elseif iter == 3, p(1) = c3;
% % %             elseif iter == 4, p(1) = c4;
% % %             end
% % %             %end
% % %             set(h1,'Position',p);
% % %             % Set title
% % %             if ftype2 == 2, title('4D Flow-STB VA'); 
% % %             elseif ftype2 == 3, title('4D Flow-CFD VA');
% % %             elseif ftype2 == 4, title('4D Flow-STB');
% % %             else, title('4D Flow-CFD'); 
% % %             end
% % %             % Save info in array
% % %             ba_results{q,1}(:,iter) = [meanDiff;CR(1);CR(2)];
% % % 
% % % 
% % %             % Save the disparity error in an array
% % %             dist_disparity_error{iter,1} = dist_disparity;
% % %             dist_disparity_avgerror(iter,1) = mean(dist_disparity);        
% % % 
% % %             % Increment the iteration counter
% % %             iter = iter + 1;    
% % %         end
% % %     end
% % %     % Get the cycle point print info
% % %     if q == 1, cycle_pt = 'TAWSS_';
% % %     elseif q == 2, cycle_pt = 'OSI_';
% % %     else, cycle_pt = 'RRT_'; 
% % %     end   
% % %     % Print the Bland-Altman plots
% % %     print(300+q,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-4DFlow-comp.pdf'])
% % % end


%%% BLAND-ALTMAN OF TAWSS, OSI, RRT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This analysis is done only comparing STB/CFD voxel vs full resolution
%%% PLOT BLAND-ALTMAN COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plotting is done at peak systole, peak diastole, and mid diastole in 
% order to compare all modalities
%mean_color = [135,206,250]/255;
mean_color = [65,105,225]/255;
mean_color_2 = [65,161,225]/255;
% Iterate through the three TAWSS, OSI, RRT
for q = 1:1:3
    % Set the parameters based on the q
    if q == 1
        u_all = tawss_all;
        x_all = tawss_all_X;
        y_all = tawss_all_Y;
        z_all = tawss_all_Z;
    elseif q == 2
        u_all = osi_all;
        x_all = osi_all_X;
        y_all = osi_all_Y;
        z_all = osi_all_Z;
    else
        u_all = rrt_all;
        x_all = rrt_all_X;
        y_all = rrt_all_Y;
        z_all = rrt_all_Z;
    end
    % Iterate through all file types
    iter = 1; % Set the iteration counter
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
    % Subplot location settings
    dc = 0.1847; c2 = 0.13+dc; c3 = 0.13+2*dc; c4 = 0.13+3*dc;
    dr = 0.25; r2 = 0.5025; %0.7093-dr; r3 = 0.7093-2*dr;
    %umin = -0.25; umax = 0.25; % Max and min for u-vel x-axis settings
    %wmin = -0.15; wmax = 0.15; % Max and min for u-vel x-axis settings
    %ymin = -0.25; ymax = 0.25; % Max and min for y-axis settings
    if strcmp(cadat.FILE.anynum,'007')
        if q == 1
            scale_fact = 10;
            ymin = -30; ymax = 12; % Max and min for y-axis settings
            dx = 12;
            x_axis_label = 'mean (dynes/cm^{2})';
            y_axis_label = 'difference (dynes/cm^{2})';
        elseif q == 2
            scale_fact = 1;
            ymin = -0.4; ymax = 0.405; % Max and min for y-axis settings
            dx = 0.2;
            x_axis_label = 'mean';
            y_axis_label = 'difference';
        else
            scale_fact = 0.1;
            ymin = -15; ymax = 15; % Max and min for y-axis settings
            dx = 5;
            x_axis_label = 'mean (dynes/cm^{2})^{-1}';
            y_axis_label = 'difference (dynes/cm^{2})^{-1}';
        end
        xmax = [24,0.4,10];
    else
        if q == 1
            scale_fact = 10;
            ymin = -30; ymax = 12; % Max and min for y-axis settings
            dx = 12;
            x_axis_label = 'mean (dynes/cm^{2})';
            y_axis_label = 'difference (dynes/cm^{2})';
        elseif q == 2
            scale_fact = 1;
            ymin = -0.4; ymax = 0.405; % Max and min for y-axis settings
            dx = 0.2;
            x_axis_label = 'mean';
            y_axis_label = 'difference';
        else
            scale_fact = 0.1;
            ymin = -15; ymax = 15; % Max and min for y-axis settings
            dx = 4;
            x_axis_label = 'mean (dynes/cm^{2})^{-1}';
            y_axis_label = 'difference (dynes/cm^{2})^{-1}';
        end
        xmax = [24,0.4,8];
    end
    
    %umin = -0.4; umax = 0.4; % Max and min for u-vel x-axis settings
    for ftype1 = 1
        % Set the data for the first file type
        x1 = x_all(ftype1,1); x1 = x1{1,1};
        y1 = y_all(ftype1,1); y1 = y1{1,1};
        z1 = z_all(ftype1,1); z1 = z1{1,1};
        u1 = u_all(ftype1,1); u1 = u1{1,1};
        
        % Set the aneurysm mask points for the first time point
        if ftype1 == 1
            % Set the aneurysm mask for the first dataset
            any_mask_1 = mri_any_mask;
            % Set the velocity mask information for the current dataset
            xmask_1 = mri_maskX;
            ymask_1 = mri_maskY;
            zmask_1 = mri_maskZ;
        elseif ftype1 == 2
            any_mask_1 = stb_any_mask;
            % Set the velocity mask information for the current dataset
            xmask_1 = stb_maskX;
            ymask_1 = stb_maskY;
            zmask_1 = stb_maskZ;
        else
            any_mask_1 = cfd_any_mask;
            % Set the velocity mask information for the current dataset
            xmask_1 = cfd_maskX;
            ymask_1 = cfd_maskY;
            zmask_1 = cfd_maskZ;
        end
        
        % Iterate through all other modalities to compare
        for ftype2 = [2,4,3,5]
            % Set the minimum threshold
            if ftype1 == 2 || ftype2 == 2 || ftype1 == 3 || ftype2 == 3
                grid_half = 0.67; %0.15;
                dist_thresh = sqrt(3*grid_half^2);
            else
                grid_half = 0.67;
                dist_thresh = sqrt(3*grid_half^2);
            end
            
            % Set the aneurysm mask points for the first time point
            if ftype2 == 4 || ftype2 == 2
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
            
            % Set the primary array as the shorter one, secondary as longer
            % Set the aneurysm mask to be used based on the primary dataset
            if length(x1)>length(x2)
                xp = x2; yp = y2; zp = z2; up = u2; 
                xs = x1; ys = y1; zs = z1; us = u1;
                any_mask = any_mask_2;
                xmask = xmask_2; ymask = ymask_2; zmask = zmask_2;
                mult_factor = -1;
            else
                xs = x2; ys = y2; zs = z2; us = u2;
                xp = x1; yp = y1; zp = z1; up = u1;
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
                in_any_mask = min_pt_dist < grid_half; % Determine if the point is in the aneurysm mask
                
                % If the minimum distance is below the threshold and the 
                % point is within the defined ROI, a match is accepted, add
                % the point to the output arrays
                if (min_dist < dist_thresh) && all_in_range && in_any_mask
                    if q == 3
                        if up(pn) < 800 && us(sn) < 800
                            dist_disparity(ct,1) = min_dist;
                            ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
                            ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
                            ct = ct + 1;
                        end
                    else
                        dist_disparity(ct,1) = min_dist;
                        ba_vec1(ct,:) = [cx,cy,cz,up(pn)];
                        ba_vec2(ct,:) = [xs(sn),ys(sn),zs(sn),us(sn)];
                        ct = ct + 1;
                    end
                end
            end
            
            % Set index skip for plotting
            iskip = 1; %if iter < 3, iskip = 1;
            %else iskip = 10; end
            
            % Compute the B-A variables and  plot them
            pbrat = 1; % Size ratio of plots
            msize = 1.5; % Set the marker size
            line_size = 2; % Set the line width size
            font_size = 10; % Set the font size

            %%% plot B-A analysis
            ba_vec1(:,4) = scale_fact*ba_vec1(:,4);
            ba_vec2(:,4) = scale_fact*ba_vec2(:,4);
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
            figure(310+q); h1 = subplot(1,4,iter);
            hold off; plot(ba_mean(1:iskip:end),ba_diff(1:iskip:end),'ok','MarkerSize',msize,'MarkerFaceColor','k'); hold on;
            plot([-20 0 100], [CR(1) CR(1) CR(1)],'r-','LineWidth',line_size); %%%plot the upper CR
            plot([-20 0 100], [CR(2) CR(2) CR(2)],'r-','LineWidth',line_size); %%%plot the lower CR
            %plot([-20 20], [-20 20].*linFit(1)+linFit(2),'--','Color',mean_color,'LineWidth',line_size); %%%plot the linear fit
            plot([-20 0 100], [meanDiff meanDiff meanDiff],'-','Color',mean_color_2,'LineWidth',line_size); %%%plot the linear fit
            %text(0.7*xmax(q),1.15*CR(1),['+CI = ',num2str(CR(1),'%.2f')]);
            %text(0.7*xmax(q),1.15*CR(2),['-CI = ',num2str(CR(2),'%.2f')]);
            %text(0.7*xmax(q),max((1.2*meanDiff),meanDiff+0.03),['MD = ',num2str(meanDiff,'%.2f')]);
            axis([0 xmax(q) ymin ymax]); set(gca,'FontSize',font_size); set(gcf,'Color',[1 1 1]); pbaspect([1 pbrat 1]);
            set(gca,'XTick',0:dx:xmax(q)); grid('on');
            if iter > 1
                set(gca,'YTickLabel',{});
            else
                ylabel(y_axis_label)
            end
            % Set the axis information
            xlabel(x_axis_label)


            p = get(h1,'Position');
            % Set the subplot location (row)
            if iter == 2, p(1) = c2;
            elseif iter == 3, p(1) = c3;
            elseif iter == 4, p(1) = c4;
            end
            %end
            set(h1,'Position',p);
            % Set title
            if ftype2 == 2, title('4D Flow-STB VA'); 
            elseif ftype2 == 3, title('4D Flow-CFD VA');
            elseif ftype2 == 4, title('4D Flow-STB');
            else, title('4D Flow-CFD'); 
            end
            % Save info in array
            ba_results{q,1}(:,iter) = [meanDiff;CR(1);CR(2)];
            
            
            % Save the disparity error in an array
            dist_disparity_error{iter,1} = dist_disparity;
            dist_disparity_avgerror(iter,1) = mean(dist_disparity);        
            
            % Increment the iteration counter
            iter = iter + 1;    
        end
    end
    % Get the cycle point print info
    if q == 1, cycle_pt = 'TAWSS_';
    elseif q == 2, cycle_pt = 'OSI_';
    else, cycle_pt = 'RRT_'; 
    end   
    % Print the Bland-Altman plots
    print(310+q,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-4DFlow-anymask-comp.pdf'])
    savefig(310+q,[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_BA_',cycle_pt,'all-4DFlow-anymask-comp.fig'])
end



if PLOT_OSI_SURFACE
    %%% Interpolate the color scheme for the WSS contours %%%%%%%%%%%%%%%%%%%%%
    % For each vertex interpolate the WSS value based on the wall shear stress
    % values
    min_pts = 3;

    % Iterate through all modalities and resolutions
    for ftype = 1:1:5
        if ftype == 1
            % Set the current x, y, z, and wss vectors
            curr_x = osi_all_X{ftype,1};
            curr_y = osi_all_Y{ftype,1};
            curr_z = osi_all_Z{ftype,1};
            curr_wss = osi_all{ftype,1};
            curr_vert = vert_mri;
            min_dist = 1.25;
        elseif ftype == 2
            % Set the current x, y, z, and wss vectors
            curr_x = osi_all_X{4,1};
            curr_y = osi_all_Y{4,1};
            curr_z = osi_all_Z{4,1};
            curr_wss = osi_all{4,1};
            curr_vert = vert_stb;
            min_dist = 0.45;
        elseif ftype == 3
            % Set the current x, y, z, and wss vectors
            curr_x = osi_all_X{5,1};
            curr_y = osi_all_Y{5,1};
            curr_z = osi_all_Z{5,1};
            curr_wss = osi_all{5,1};
            curr_vert = vert_cfd;
            min_dist = 0.45;
        elseif ftype == 4
            % Set the current x, y, z, and wss vectors
            curr_x = osi_all_X{2,1};
            curr_y = osi_all_Y{2,1};
            curr_z = osi_all_Z{2,1};
            curr_wss = osi_all{2,1};
            curr_vert = vert_stb;
            min_dist = 1.25;
        else
            % Set the current x, y, z, and wss vectors
            curr_x = osi_all_X{3,1};
            curr_y = osi_all_Y{3,1};
            curr_z = osi_all_Z{3,1};
            curr_wss = osi_all{3,1};
            curr_vert = vert_cfd;
            min_dist = 1.25;
        end

        % Initialize the output color matrix and interpolate the color
        curr_osi_color = zeros(size(curr_vert,1),1);
        for q = 1:1:size(curr_vert,1)
            % Get the current vertex 
            cx = curr_vert(q,1);
            cy = curr_vert(q,2);
            cz = curr_vert(q,3);

            % Get the closest WSS points to the current vertex
            pt_dists = sqrt((curr_x-cx).^2 + (curr_y-cy).^2 + (curr_z-cz).^2);
            [srt_dists,srt_pts] = sort(pt_dists,'ascend');
            srt_wss = curr_wss(srt_pts);

            % Identify the points which are to be included in the interpolated
            % color calculation
            curr_dist = min_dist; % Set the current distance
            num_pts = sum(srt_dists <= curr_dist); % Compute the number of points in the minimum distance
            while num_pts < min_pts % Iterate until the minimum number of points are included
                curr_dist = curr_dist + 0.1; % Increase the distance 
                num_pts = sum(srt_dists <= curr_dist); % Recompute the number of points included
            end
            % Collect the indices to be included in the calculation
            inc_inds = (srt_dists <= curr_dist); % Gather the indices of points to include
            inc_wss = srt_wss(inc_inds == 1); % Gather the WSS of points to include
            inc_dists = srt_dists(inc_inds == 1); % Gather the distances of points to include
            inc_wgts = 1./(inc_dists.^2); % Compute the weights
            inc_wgts = inc_wgts/sum(inc_wgts); % Normalize the weights
            interp_osi = sum(inc_wgts.*inc_wss);

            % Save the interpolated WSS in the wss color vector
            curr_osi_color(q) = interp_osi;
        end
        % Set the current color to the appropriate array
        if ftype == 1
            mri_wss_color = curr_osi_color;
        elseif ftype == 2
            stb_fres_wss_color = curr_osi_color; 
        elseif ftype == 3
            cfd_fres_wss_color = curr_osi_color;
        elseif ftype == 4
            stb_wss_color = curr_osi_color; 
        else
            cfd_wss_color = curr_osi_color;
        end
    end

    %%% Plot the WSS contours %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = load('wss_cmap.mat');
    wss_cmap = fp.cmap;
    max_color = 0.5;
    fsize = 28;

    % Plot the surface contours
    figure(410); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_mri,'Faces',face_mri,'FaceVertexCData',mri_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('MRI'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;

    figure(411); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_stb,'Faces',face_stb,'FaceVertexCData',stb_fres_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Full Res'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;

    figure(412); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_cfd,'Faces',face_cfd,'FaceVertexCData',cfd_fres_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Full Res'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;

    figure(413); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_stb,'Faces',face_stb,'FaceVertexCData',stb_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('STB - Voxel Avg'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;

    figure(414); hold off;
    %[faces,verts] = isosurface(xm,ym,zm,vv,0,cfd_wssmag);
    p2 = patch('Vertices',vert_cfd,'Faces',face_cfd,'FaceVertexCData',cfd_wss_color,...
        'FaceColor','interp','EdgeColor','none'); colormap(gca,wss_cmap); caxis([0 max_color])
    p2.FaceAlpha = 0.5;
    set(gca,'FontSize',fsize); set(gcf,'Color',[1 1 1]); title('CFD - Voxel Avg'); grid(gca,'on')
    axis(axis_range); view(view_angle); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); colorbar;


    print(410,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_mri_osi_contour.pdf'])
    print(411,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_fres_osi_contour.pdf'])
    print(412,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_fres_osi_contour.pdf'])
    print(413,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_stb_osi_contour.pdf'])
    print(414,'-dpdf',[cadat.DIR.saveplots,'ANY-',cadat.FILE.anynum,'_cfd_osi_contour.pdf'])
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



end



