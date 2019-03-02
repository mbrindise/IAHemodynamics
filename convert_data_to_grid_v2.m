clear params
%%% SELECT OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select aneurysm number and file type
params.file.anynum = '007';
params.file.filetype = 'STB';
params.system.cpu = 'server'; % Can be pc or server
filetype_dt = 1/2000;%0.0015;
% Data import options
params.options.STARTTIME = 1;
params.options.LASTTIME = 898;
params.options.numtime = 3;
% STL Information
params.stl.stl_save_name = ['STL-ANY-',params.file.anynum];
% Velocity information
params.file.run_num = '1';
if strcmp(params.file.filetype,'STB')
    params.file.basename = ['any',params.file.anynum,'_stb_v',params.file.run_num,'_'];
else
    params.file.basename = ['any',params.file.anynum,'_',params.file.filetype,'_'];
end
params.file.velname = [params.file.basename,'regmask_'];
params.file.outvelname = [params.file.basename,'regmask_gridded_'];

%% DIRECTORY AND FILE NAME INFORMATION
% Set directory info based on system
if strcmp(params.system.cpu,'pc')
    params.system.file_base = 'Z:\Projects\Cerebral_Aneurysm\ANY\';
    fslash = '\';
else
    params.system.file_base = '/home/shannon/a/mbrindis/Projects/Cerebral_Aneurysm/ANY/';
    fslash = '/';
end
% Directory information
fpath = [params.system.file_base,params.file.anynum,fslash,params.file.filetype,fslash,'vel_registered_masked',fslash];
savefiles = [params.system.file_base,params.file.anynum,fslash,'registration_files',fslash];
savevel = [params.system.file_base,params.file.anynum,fslash,params.file.filetype,fslash,'vel_gridded',fslash];
savemask = [params.system.file_base,params.file.anynum,fslash,'registration_files',fslash,'ANY-',params.file.anynum,'_',params.file.filetype,'_masked_velocity_grid.mat'];

%%% CREATE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(params.file.filetype,'MRI')
    % Load an MRI velocity file
    fp = load([savevel,params.file.outvelname,num2str(1,'%05i'),'.mat']);
    xmask = unique(fp.x);
    ymask = unique(fp.y);
    zmask = unique(fp.z);
    velmask = nan(length(ymask),length(xmask),length(zmask));
    % Set the velocity mask to zero if a vector is at that location or NaN
    % if no vector is at that location (it is outside the mask)
    for q = 1:1:length(fp.x)
        % Current vector point
        cx = fp.x(q);
        cy = fp.y(q);
        cz = fp.z(q);
        % Find the mask indices it is closest to
        [~,r_ind] = min(abs(ymask-cy));
        [~,c_ind] = min(abs(xmask-cx));
        [~,z_ind] = min(abs(zmask-cz));
        % Set that mask point to 0
        velmask(r_ind,c_ind,z_ind) = 0;
    end
    % Save the mask grid
    save(savemask,'velmask','xmask','ymask','zmask')
else
    % Grid conversion options
    params.grid.dx = 0.4;
    params.grid.dy = 0.4;
    params.grid.dz = 0.4;

    % Minimum and maximum of X, Y, Z of particles
    params.grid.minX = -15;% -20;%-10; %min(params.data.X);
    params.grid.maxX = 8; %20;%8; %8; %44; %max(params.data.X);
    params.grid.minY = -6;%-15;%-5;%-8; %-5; %-16; %min(params.data.Y);
    params.grid.maxY = 24;%15;%22; %15; %30; %max(params.data.Y);
    params.grid.minZ = -12;%-7;%-8; %-7; %-14; %min(params.data.Z);
    params.grid.maxZ = 15;%7;%5;%16; %7; %12; %max(params.data.Z);

    % Form x, y, z grid based on minimum and maximum X, Y, Z of particles
    params.grid.x = (params.grid.minX-0.5*params.grid.dx):params.grid.dx:(params.grid.maxX+params.grid.dx);
    params.grid.y = (params.grid.minY-0.5*params.grid.dy):params.grid.dy:(params.grid.maxY+params.grid.dy);
    params.grid.z = (params.grid.minZ-0.5*params.grid.dz):params.grid.dz:(params.grid.maxZ+params.grid.dz);

    % Convert gird to mesh grid and to a linear version of the mesh grid
    [params.grid.X,params.grid.Y,params.grid.Z] = meshgrid(params.grid.x,params.grid.y,params.grid.z);
    params.grid.Xl = reshape(params.grid.X,[length(params.grid.x)*length(params.grid.y)*length(params.grid.z),1]);
    params.grid.Yl = reshape(params.grid.Y,[length(params.grid.x)*length(params.grid.y)*length(params.grid.z),1]);
    params.grid.Zl = reshape(params.grid.Z,[length(params.grid.x)*length(params.grid.y)*length(params.grid.z),1]);

    % Load the STL mask
    params.stl.mask = load([savefiles,params.stl.stl_save_name,'_',params.file.filetype,'_STLMask_filt.mat']);
    xm = params.stl.mask.stl_x;
    ym = params.stl.mask.stl_y;
    zm = params.stl.mask.stl_z;
    stlmask = params.stl.mask.stlmaskf;

    % Initialize the velocity mask grid
    %velmask = zeros(size(params.grid.Xl));

    % Keep only the grid points that are in the STL mask, NaN the rest
    ubase = zeros(size(params.grid.Xl));
    for q = 1:1:length(ubase)
        % Get the current x,y,z coordiantes
        xcurr = params.grid.Xl(q);
        ycurr = params.grid.Yl(q);
        zcurr = params.grid.Zl(q);

        % Find the mask point that is closest to the grid point
        [~,xind] = min(abs(xm - xcurr));
        [~,yind] = min(abs(ym - ycurr));
        [~,zind] = min(abs(zm - zcurr));

        % Determine if the current grid point is in the mask
        in_mask = stlmask(yind,xind,zind);

        % Ensure the current point is not larger than the stl mask parameters
        in_x = (xcurr >= min(xm)) & (xcurr <= max(xm));
        in_y = (ycurr >= min(ym)) & (ycurr <= max(ym));
        in_z = (zcurr >= min(zm)) & (zcurr <= max(zm));
        if (in_x == 0) || (in_y == 0) || (in_z == 0)
            in_mask = 0;
        end

        if in_mask == 0
            ubase(q) = NaN;
        end
    end
    % Convert ubase into the velocity mask grid
    velmask = reshape(ubase,[length(params.grid.y),length(params.grid.x),length(params.grid.z)]);
    xmask = params.grid.x;
    ymask = params.grid.y;
    zmask = params.grid.z;
    save(savemask,'velmask','xmask','ymask','zmask')

    % Limit the grid points to only the ones in the mask
    inmask_inds = find(~isnan(ubase));
    params.grid.xm = params.grid.Xl(inmask_inds);
    params.grid.ym = params.grid.Yl(inmask_inds);
    params.grid.zm = params.grid.Zl(inmask_inds);

    %%% READ IN MAT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_starts = params.options.STARTTIME:params.options.numtime:params.options.LASTTIME-params.options.numtime;
    dt = params.options.numtime * filetype_dt;
    for curr_t = 1:1:length(time_starts) 
        tic
        curr_act_time = time_starts(curr_t);
        fprintf('\nNow processing %i ...',curr_act_time)
        % Load files for current time step
        curr_t_x = [];
        curr_t_y = [];
        curr_t_z = [];
        curr_t_u = [];
        curr_t_v = [];
        curr_t_w = [];
        t = (curr_t-1)*dt;
        for q = 0:1:params.options.numtime-1
            fp = load([fpath,params.file.velname,num2str(time_starts(curr_t)+q,'%05i'),'.mat']);
            new_x = fp.x; % X is in mm
            new_y = fp.y; % Y is in mm
            new_z = fp.z; % Z is in mm
            new_u = fp.u; % U is in m/s
            new_v = fp.v; % V is in m/s
            new_w = fp.w; % W is in m/s
            % Ensure all vectors are a column vector
            if size(new_x,2)>1, new_x = new_x'; end
            if size(new_y,2)>1, new_y = new_y'; end
            if size(new_z,2)>1, new_z = new_z'; end
            if size(new_u,2)>1, new_u = new_u'; end
            if size(new_v,2)>1, new_v = new_v'; end
            if size(new_w,2)>1, new_w = new_w'; end
            % add the new vectors to the current one
            curr_t_x = [curr_t_x;new_x];
            curr_t_y = [curr_t_y;new_y];
            curr_t_z = [curr_t_z;new_z];
            curr_t_u = [curr_t_u;new_u];
            curr_t_v = [curr_t_v;new_v];
            curr_t_w = [curr_t_w;new_w];
        end
        fprintf(' imported ')

        %%% CONVERT TO STRUCTURED GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each grid point, identify particles within an adaptive radius.
        % Need minimum of 3 particles. Base radius is half the distance between
        % grid points.
        % Initialize 3D structure for time step
        cu = zeros(size(params.grid.xm));
        cv = zeros(size(params.grid.xm));
        cw = zeros(size(params.grid.xm));
        start_rad = 0.5*params.grid.dx;
        if strcmp(params.file.filetype,'CFD')
            max_rad = sqrt(((0.25*params.grid.dx)^2+(0.25*params.grid.dy)^2+(0.25*params.grid.dz)^2));
        else
            max_rad = 1e5*start_rad;
        end
        for p = 1:1:length(params.grid.xm)
            % Obtain x, y, z of current grid points
            cx = params.grid.xm(p);
            cy = params.grid.ym(p);
            cz = params.grid.zm(p);

            % Compute radial distance for all particles in the field
            rad_dist = sqrt((curr_t_x - cx).^2 + (curr_t_y - cy).^2 + (curr_t_z - cz).^2);

            % Find all particles within a radius of 1/2 distance of the grid
            % displacement
            %min_rad_dist(p,1) = min(rad_dist);
            % Find all points within a given radius
            pts_in_rad = (rad_dist <= start_rad); % All points within starting radius
            if sum(pts_in_rad) < 3 || (curr_rad > max_rad)
                % If not enough points in radius, need to adapt radius
                curr_rad = min(rad_dist);
                while sum(pts_in_rad) < 3
                    curr_rad = curr_rad + 0.01;
                    pts_in_rad = (rad_dist <= curr_rad); % All points within starting radius
                end
            end
            if curr_rad <= max_rad
                % Get velocity and distance for points in radius
                in_rad_u = curr_t_u(pts_in_rad == 1);
                in_rad_v = curr_t_v(pts_in_rad == 1);
                in_rad_w = curr_t_w(pts_in_rad == 1);
                in_rad_dists = rad_dist(pts_in_rad == 1);

                % Compute weight based on square of distance
                dist_wgt = 1./(in_rad_dists.^2);
                dist_wgt = dist_wgt/sum(dist_wgt);

                % Save weighted velocity for point
                cu(p) = sum(in_rad_u.*dist_wgt);
                cv(p) = sum(in_rad_v.*dist_wgt);
                cw(p) = sum(in_rad_w.*dist_wgt);
                rad_save(p) = curr_rad;
            else
                cu(p) = nan;
                cv(p) = nan;
                cw(p) = nan;
            end
        end
        fprintf('- gridded')
        % Keep only indices whose velocities are not NaN values
        if curr_t == 1
            kp_inds = ~isnan(cu);
        end
        
        % Get the output vectors
        x = params.grid.xm(kp_inds);
        y = params.grid.ym(kp_inds);
        z = params.grid.zm(kp_inds);
        u = cu(kp_inds);
        v = cv(kp_inds);
        w = cw(kp_inds);
        
        % Update the velocity mask 
        if curr_t == 1
            ubase = zeros(size(params.grid.Xl));
            for q = 1:1:length(ubase)
                % Get the current x,y,z coordiantes
                xcurr = params.grid.Xl(q);
                ycurr = params.grid.Yl(q);
                zcurr = params.grid.Zl(q);

                % Find the mask point that is closest to the grid point
                all_pt_dist = sqrt((x-xcurr).^2 + (y-ycurr).^2 + (z-zcurr).^2);
                pt_dist = min(all_pt_dist);
                close_pt = pt_dist < 0.01;

                % Determine if the current grid point is in the mask
                % Find the mask point that is closest to the grid point
                [~,xind] = min(abs(xm - xcurr));
                [~,yind] = min(abs(ym - ycurr));
                [~,zind] = min(abs(zm - zcurr));
                in_mask = stlmask(yind,xind,zind);

                % Ensure the current point is not larger than the stl mask parameters
                in_x = (xcurr >= min(xm)) & (xcurr <= max(xm));
                in_y = (ycurr >= min(ym)) & (ycurr <= max(ym));
                in_z = (zcurr >= min(zm)) & (zcurr <= max(zm));
                if (in_x == 0) || (in_y == 0) || (in_z == 0)
                    in_mask = 0;
                end
                if close_pt == 0
                    in_mask = 0;
                end          
                if in_mask == 0
                    ubase(q) = NaN;
                end
            end
            % Convert ubase into the velocity mask grid
            velmask = reshape(ubase,[length(params.grid.y),length(params.grid.x),length(params.grid.z)]);
            xmask = params.grid.x;
            ymask = params.grid.y;
            zmask = params.grid.z;
            save(savemask,'velmask','xmask','ymask','zmask')
        end

        % Save the grid file in MAT and DAT formats
        % Save mat file
        save([savevel,params.file.outvelname,num2str(curr_act_time,'%05i'),'.mat'],'x','y','z','u','v','w','t');

        % Save dat file
        if curr_act_time <= 30
            % Initialize output file
            fid = fopen([savevel,params.file.outvelname,num2str(curr_act_time,'%05i'),'.dat'],'w');
            fprintf(fid,['TITLE = "',params.file.outvelname,num2str(curr_act_time,'%05i'),'"']);
            fprintf(fid,'\r\nVARIABLES = "X", "Y", "Z", "V-X", "V-Y", "V-Z"');

            % Iterate through all points and write only those that are non-zero
            for zz = 1:1:length(x)
                cx = x(zz);
                cy = y(zz);
                cz = z(zz);
                cu = u(zz);
                cv = v(zz);
                cw = w(zz);
                fprintf(fid,'\r\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
            end
            % Close dat file
            fclose(fid);
        end
        tval = toc;
        tmin = floor(tval/60);
        tsec = round(tval - tmin*60);
        fprintf(' - completed successfully in %01i:%02i',tmin,tsec)
    end
    fprintf('\n\n')
end




