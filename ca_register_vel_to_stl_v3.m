%%% CEREBRAL ANEURYSM VELOCITY REGISTRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program registers the input velocity fields using the STL file. All
% inputs need to be logged in the cadat structure. This structure can be
% formed using the create_cadat_structure script. All needed information is
% indicated in this script. The registration process should only need to be
% run 1 time.
%
% This program requires user inputs to determine if certain steps are
% needed.
%
% The output of this program is:
% (1) The registered and masked x, y coordinates
% (2) The masked velocity fields
%
% These two outputs can be directly inputted into the convert_data_to_grid
% program which converts CFD and STB data to a grid. For tomo PIV data,
% this step is not needed as a grid is already formed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cadat] = ca_register_vel_to_stl_v3(cadat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARE STL AND VELOCITY FIELD DATA, FLIP DATA AS NECESSARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in data velocity field, could be PIV, CFD, or MRI data %%%%%%%%%%%
% STL/velocity flipping to match fields will be decided using first time step
for vnum = cadat.TIME.ts:1:cadat.TIME.te
    fprintf('---------- Now processing file %i  ----------',vnum)
    fprintf('\nImporting velocity file %i ...',vnum)
    if strcmp(cadat.FILE.filetype,'TOMO')
        fid = fopen([cadat.DIR.velfiles,cadat.FILE.basename,num2str(vnum,'%05i'),'.dat'],'rt');
        A = textscan(fid,'%f %f %f %f %f %f','HeaderLines',3);
        cadat.DATA.x = A{1,1};
        cadat.DATA.y = A{1,2};
        cadat.DATA.z = A{1,3};
        cadat.DATA.u = A{1,4};
        cadat.DATA.v = A{1,5};
        cadat.DATA.w = A{1,6};

        % Filter velocity fields by removing all zero velocities
        num_vel_pts = size(cadat.DATA.x,1);
        ct = 1;

        for zz = 1:1:num_vel_pts
            if cadat.DATA.u(zz) ~= 0 && cadat.DATA.v(zz) ~= 0 && cadat.DATA.w(zz) ~= 0
                cadat.DATA.xf(ct,1) = cadat.DATA.x(zz);
                cadat.DATA.yf(ct,1) = cadat.DATA.y(zz);
                cadat.DATA.zf(ct,1) = cadat.DATA.z(zz);
                cadat.DATA.uf(ct,1) = cadat.DATA.u(zz);
                cadat.DATA.vf(ct,1) = cadat.DATA.v(zz);
                cadat.DATA.wf(ct,1) = cadat.DATA.w(zz);
                ct = ct + 1;
            end
        end
        
        % Create mesh grid using dx and dy from velocity field
        unq_x = unique(cadat.DATA.xf);
        unq_y = unique(cadat.DATA.yf);
        unq_z = unique(cadat.DATA.zf);
        cadat.DATA.dx = mean(unq_x(2:end)-unq_x(1:end-1));
        cadat.DATA.dy = mean(unq_y(2:end)-unq_y(1:end-1));
        cadat.DATA.dz = mean(unq_z(2:end)-unq_z(1:end-1));
    else
        % If STB is used, load the mat file, make it the filtered velocity
        % because there is no need to remove zeros
        A = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(vnum,'%06i'),'.mat']);
        cadat.DATA.x = A.x;
        cadat.DATA.y = A.y;
        cadat.DATA.z = A.z;
        cadat.DATA.u = A.u;
        cadat.DATA.v = A.v;
        cadat.DATA.w = A.w;
        cadat.DATA.xf = A.x;
        cadat.DATA.yf = A.y;
        cadat.DATA.zf = A.z;
        cadat.DATA.uf = A.u;
        cadat.DATA.vf = A.v;
        cadat.DATA.wf = A.w;
        
        % For STB, dx, dy, dz needs to be specified already - no
        % computation necessary
        
        % If the input files are not in mm, needs to be converted to mm
        if cadat.DATA.inMM == 0
            cadat.DATA.x = cadat.DATA.x * 1000; % Convert to mm
            cadat.DATA.y = cadat.DATA.y * 1000;
            cadat.DATA.z = cadat.DATA.z * 1000;
            cadat.DATA.xf = cadat.DATA.x; % Convert to mm
            cadat.DATA.yf = cadat.DATA.y;
            cadat.DATA.zf = cadat.DATA.z;
        end
        
        % Keep only the data that has non-zero velocities
        vel_kp = (abs(cadat.DATA.u) + abs(cadat.DATA.v) + abs(cadat.DATA.w)) > 0;
        cadat.DATA.x = cadat.DATA.x(vel_kp == 1);
        cadat.DATA.y = cadat.DATA.y(vel_kp == 1);
        cadat.DATA.z = cadat.DATA.z(vel_kp == 1);
        cadat.DATA.u = cadat.DATA.u(vel_kp == 1);
        cadat.DATA.v = cadat.DATA.v(vel_kp == 1);
        cadat.DATA.w = cadat.DATA.w(vel_kp == 1);
        cadat.DATA.xf = cadat.DATA.x;
        cadat.DATA.yf = cadat.DATA.y;
        cadat.DATA.zf = cadat.DATA.z;
        cadat.DATA.uf = cadat.DATA.u;
        cadat.DATA.vf = cadat.DATA.v;
        cadat.DATA.wf = cadat.DATA.w;
    end

    % Get x and y line coordinates
    cadat.DATA.Xl = min(cadat.DATA.x):cadat.DATA.dx:(max(cadat.DATA.x) + 0.25*cadat.DATA.dx);
    cadat.DATA.Yl = min(cadat.DATA.y):cadat.DATA.dy:(max(cadat.DATA.y) + 0.25*cadat.DATA.dy);
    cadat.DATA.Zl = min(cadat.DATA.z):cadat.DATA.dz:(max(cadat.DATA.z) + 0.25*cadat.DATA.dz);
    cadat.DATA.Xlc = cadat.DATA.Xl - 0.5*(max(cadat.DATA.Xl) + min(cadat.DATA.Xl));
    cadat.DATA.Ylc = cadat.DATA.Yl - 0.5*(max(cadat.DATA.Yl) + min(cadat.DATA.Yl));
    cadat.DATA.Zlc = cadat.DATA.Zl - 0.5*(max(cadat.DATA.Zl) + min(cadat.DATA.Zl));
    
    % Save the center values for the first time step - this must be the
    % time step used for the shifting
    if vnum == cadat.TIME.ts
        cadat.DATA.centvalX = 0.5*(max(cadat.DATA.Xl) + min(cadat.DATA.Xl));
        cadat.DATA.centvalY = 0.5*(max(cadat.DATA.Yl) + min(cadat.DATA.Yl));
        cadat.DATA.centvalZ = 0.5*(max(cadat.DATA.Zl) + min(cadat.DATA.Zl));
    end

    % Create grid
    [cadat.DATA.X,cadat.DATA.Y,cadat.DATA.Z] = meshgrid(cadat.DATA.Xlc,cadat.DATA.Ylc,cadat.DATA.Zlc);

    % Reshape velocity into grid format
    cadat.DATA.U = zeros(size(cadat.DATA.X));
    cadat.DATA.V = zeros(size(cadat.DATA.X));
    cadat.DATA.W = zeros(size(cadat.DATA.X));

    for zz = 1:1:length(cadat.DATA.xf)
        cx = cadat.DATA.xf(zz);
        cy = cadat.DATA.yf(zz);
        cz = cadat.DATA.zf(zz);
        cu = cadat.DATA.uf(zz);
        cv = cadat.DATA.vf(zz);
        cw = cadat.DATA.wf(zz);

        [~,zind] = min(abs(cadat.DATA.Zl - cz));
        [~,cind] = min(abs(cadat.DATA.Xl - cx));
        [~,rind] = min(abs(cadat.DATA.Yl - cy));

        cadat.DATA.U(rind,cind,zind) = cu;
        cadat.DATA.V(rind,cind,zind) = cv;
        cadat.DATA.W(rind,cind,zind) = cw;

    end
    
    fprintf('completed successfully')

    %%% LOAD IN STL POINT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This step only needs to be performed 1 single time. Thus if it is
    % already created in previous time step/run of code, the adjusted stl
    % file only needs to be loaded.
    % Each modality should also use the SAME adjusted STL! So this step
    % only needs to happen 1 time for 1 modality!
    if (cadat.OPTIONS.ADJSTL == 1) && (vnum == cadat.TIME.ts) && (cadat.OPTIONS.TESTRESOLUTION == 0)
        fprintf('\nAdjusting STL file...\n')
        % Read STL file
        [cadat.STL.v, cadat.STL.f, cadat.STL.n, cadat.STL.c, ~ ] = stlread([cadat.DIR.stlfiles,cadat.STL.fp_stl_name,'.stl']);
        cadat.STL.x = cadat.STL.v(:,1);
        cadat.STL.y = cadat.STL.v(:,2);
        cadat.STL.z = cadat.STL.v(:,3);
        
        % Need to determine if any coordinates need to be switched
        %close all;
        figure(20); hold off; scatter3(cadat.STL.x,cadat.STL.y,cadat.STL.z);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        flip_xy = input('Flip X and Y coordinates? [1-Yes, 0-No]: ');
        if flip_xy
            tempvar = cadat.STL.x;
            cadat.STL.x = cadat.STL.y;
            cadat.STL.y = tempvar;
        end
        figure(20); hold off; scatter3(cadat.STL.x,cadat.STL.y,cadat.STL.z);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        flip_xz = input('Flip X and Z coordinates? [1-Yes, 0-No]: ');
        if flip_xz
            tempvar = cadat.STL.x;
            cadat.STL.x = cadat.STL.z;
            cadat.STL.z = tempvar;
        end
        figure(20); hold off; scatter3(cadat.STL.x,cadat.STL.y,cadat.STL.z);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        flip_yz = input('Flip Y and Z coordinates? [1-Yes, 0-No]: ');
        if flip_yz
            tempvar = cadat.STL.y;
            cadat.STL.y = cadat.STL.z;
            cadat.STL.z = tempvar;
        end
        
        % Zero the x, y, z of the STL field. Flip field if necessary
        mid_x = 0.5*(max(cadat.STL.x) + min(cadat.STL.x));
        mid_y = 0.5*(max(cadat.STL.y) + min(cadat.STL.y));
        mid_z = 0.5*(max(cadat.STL.z) + min(cadat.STL.z));
        cadat.STL.x = cadat.STL.x - mid_x;
        cadat.STL.y = cadat.STL.y - mid_y;
        cadat.STL.z = cadat.STL.z - mid_z;

        % Identify Z-planes in velocity field
        cadat.DATA.zplanes = unique(cadat.DATA.zf);
        cadat.DATA.zmax = max(cadat.DATA.zplanes);
        cadat.DATA.zmin = min(cadat.DATA.zplanes);

        % Identify Z-planes in STL file
        cadat.STL.zplanes = unique(cadat.STL.z);

        % Consolodate STL Z-planes to a chosen STL dz
        new_stl_z = min(cadat.STL.zplanes):cadat.STL.dz:max(cadat.STL.zplanes);
        if new_stl_z(end) < max(cadat.STL.zplanes) % Ensure full range met
            new_max = new_stl_z(end) + cadat.DATA.dz;
            new_stl_z = [new_stl_z,new_max];
        end
        % Place each STL point on a gridded z-plane
        for zz = 1:1:length(cadat.STL.z)
            min_z_diff = abs(new_stl_z - cadat.STL.z(zz));
            [~,best_z] = min(min_z_diff);
            cadat.STL.zadj(zz) = new_stl_z(best_z);
        end
        cadat.STL.zplanes_adj = unique(cadat.STL.zadj);

        % Remove redunant points after consolidation
        next_ind = 1;
        prev_mark = 0;
        for zz = 1:1:length(cadat.STL.zplanes_adj)
            pt2kp = [];
            curr_z = cadat.STL.zplanes_adj(zz);
            curr_z_inds = find(cadat.STL.zadj == curr_z);
            curr_x_pts = cadat.STL.x(curr_z_inds);
            curr_y_pts = cadat.STL.y(curr_z_inds);

            for z2 = 1:1:length(curr_z_inds)
                cpt_x = cadat.STL.x(curr_z_inds(z2));
                cpt_y = cadat.STL.y(curr_z_inds(z2));
                x_pts = find(curr_x_pts == cpt_x);
                y_pts = find(curr_y_pts == cpt_y);
                same_pts = intersect(x_pts,y_pts);

                pt2kp = [pt2kp,same_pts(1)];
                pt2kp = unique(pt2kp);

            end
            stl_rm_x(next_ind:next_ind+length(pt2kp)-1) = cadat.STL.x(curr_z_inds(pt2kp));
            stl_rm_y(next_ind:next_ind+length(pt2kp)-1) = cadat.STL.y(curr_z_inds(pt2kp));
            stl_rm_z(next_ind:next_ind+length(pt2kp)-1) = cadat.STL.zadj(curr_z_inds(pt2kp));
            next_ind = length(stl_rm_x) + 1;
            
            % Print the percent completed for steps of 10
            perc_complete = floor(100*zz/length(cadat.STL.zplanes_adj));
            if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                fprintf(' %i%%',perc_complete)
                prev_mark = perc_complete;
            end
        end
        save([cadat.DIR.savefiles,cadat.STL.adjstlptsfile,'.mat'],'stl_rm_x','stl_rm_y','stl_rm_z')
        fprintf(' completed successfully')
    else
        load([cadat.DIR.savefiles,cadat.STL.adjstlptsfile,'.mat'])
    end
    %%% DETERMINE FLIPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine if the velocity fields need to be flipped, or if the
    % STL files need to be flipped. This step requires user inputs. If the
    % data is to be flipped, the velocities also need to be flipped. This
    % step also only needs to be completed one time and then can be loaded
    % for the remaining time steps
    if (cadat.OPTIONS.DETERMINE_FLIPS == 1) &&  (vnum == cadat.TIME.ts) && (cadat.OPTIONS.TESTRESOLUTION == 0)
        fprintf('\nDetermining flips...\n')
        flip_field_x = 0; % Initialize flips
        flip_field_y = 0; % Initialize flips
        flip_field_z = 0; % Initialize flips
        flip_stl_x = 0; % Initialize flips
        flip_stl_y = 0; % Initialize flips
        flip_stl_z = 0; % Initialize flips
        
        % Load in the flipped STL coordinates if the flipping is not to be
        % saved - this way the new dataset type can be flipped to the STL
        if cadat.OPTIONS.save_stl_flips == 0
            fp = load([cadat.DIR.savefiles,cadat.STL.adjstlflipfile,'.mat']);
            stl_rm_x = fp.stl_rm_x;
            stl_rm_y = fp.stl_rm_y;
            stl_rm_z = fp.stl_rm_z;
        end
        
        % Determine signs of coordinate flips
        close all;
        figure(25); subplot(1,2,1); scatter3(stl_rm_x,stl_rm_y,stl_rm_z); view(2);
        set(gca,'FontSize',14); title('STL'); xlabel('X'); ylabel('Y')
        subplot(1,2,2); scatter3(cadat.DATA.x,cadat.DATA.y,cadat.DATA.z); view(2);
        set(gca,'FontSize',14); title('DATA'); xlabel('X'); ylabel('Y')
        % Run through all prompts and 
        prompt = 'Flip X (Left-Right)? 0-Neither, 1-STL only, 2-Data only, 3-Both:  ';
        flip_lr_choice = input(prompt);
        % Determine X Flip
        switch flip_lr_choice
            case 1 % only flip STL file
                stl_rm_x = -stl_rm_x;
                flip_stl_x = 1;
            case 2 % only flip data
                cadat.DATA.x = -cadat.DATA.x;
                cadat.DATA.xf = -cadat.DATA.xf;
                cadat.DATA.Xl = sort(-cadat.DATA.Xl,'ascend');
                cadat.DATA.Xlc = sort(-cadat.DATA.Xlc,'ascend');
                cadat.DATA.X = -cadat.DATA.X;
                cadat.DATA.u = -cadat.DATA.u;
                flip_field_x = 1;
            case 3 % flip both stl and data
                stl_rm_x = -stl_rm_x;
                cadat.DATA.x = -cadat.DATA.x;
                cadat.DATA.xf = -cadat.DATA.xf;
                cadat.DATA.Xl = sort(-cadat.DATA.Xl,'ascend');
                cadat.DATA.Xlc = sort(-cadat.DATA.Xlc,'ascend');
                cadat.DATA.X = -cadat.DATA.X;
                cadat.DATA.u = -cadat.DATA.u;
                flip_field_x = 1;
                flip_stl_x = 1;
        end
        % Determine Y Flip
        figure(25); subplot(1,2,1); scatter3(stl_rm_x,stl_rm_y,stl_rm_z); view(2);
        set(gca,'FontSize',14); title('STL'); xlabel('X'); ylabel('Y')
        subplot(1,2,2); scatter3(cadat.DATA.x,cadat.DATA.y,cadat.DATA.z); view(2);
        set(gca,'FontSize',14); title('DATA'); xlabel('X'); ylabel('Y')
        prompt = 'Flip Y (Up-Down)? 0-Neither, 1-STL only, 2-Data only, 3-Both:  ';
        flip_ud_choice = input(prompt);
        switch flip_ud_choice
            case 1 % only flip STL file
                stl_rm_y = -stl_rm_y;
                flip_stl_y = 1;
            case 2 % only flip data
                cadat.DATA.y = -cadat.DATA.y;
                cadat.DATA.yf = -cadat.DATA.yf;
                cadat.DATA.Yl = sort(-cadat.DATA.Yl,'ascend');
                cadat.DATA.Ylc = sort(-cadat.DATA.Ylc,'ascend');
                cadat.DATA.Y = -cadat.DATA.Y;
                cadat.DATA.v = -cadat.DATA.v;
                flip_field_y = 1;
            case 3 % flip both stl and data
                stl_rm_y = -stl_rm_y;
                cadat.DATA.y = -cadat.DATA.y;
                cadat.DATA.yf = -cadat.DATA.yf;
                cadat.DATA.Yl = sort(-cadat.DATA.Yl,'ascend');
                cadat.DATA.Ylc = sort(-cadat.DATA.Ylc,'ascend');
                cadat.DATA.Y = -cadat.DATA.Y;
                cadat.DATA.v = -cadat.DATA.v;
                flip_field_y = 1;
                flip_stl_y = 1;
        end
        % Determine Z Flip
        figure(25); subplot(1,2,1); scatter3(stl_rm_x,stl_rm_y,stl_rm_z); view(90,0);
        set(gca,'FontSize',14); title('STL'); zlabel('Z'); ylabel('Y')
        subplot(1,2,2); scatter3(cadat.DATA.x,cadat.DATA.y,cadat.DATA.z); view(90,0);
        set(gca,'FontSize',14); title('DATA'); zlabel('Z'); ylabel('Y')
        prompt = 'Flip Z (Back-Forth)? 0-Neither, 1-STL only, 2-Data only, 3-Both:  ';
        flip_z_choice = input(prompt);
        switch flip_z_choice
            case 1 % only flip STL file
                stl_rm_z = -stl_rm_z;
                flip_stl_z = 1;
            case 2 % only flip data
                cadat.DATA.z = -cadat.DATA.z;
                cadat.DATA.zf = -cadat.DATA.zf;
                cadat.DATA.Zl = sort(-cadat.DATA.Zl,'ascend');
                cadat.DATA.Zlc = sort(-cadat.DATA.Zlc,'ascend');
                cadat.DATA.Z = -cadat.DATA.Z;
                cadat.DATA.w = -cadat.DATA.w;
                flip_field_z = 1;
            case 3 % flip both stl and data
                stl_rm_z = -stl_rm_z;
                cadat.DATA.z = -cadat.DATA.z;
                cadat.DATA.zf = -cadat.DATA.zf;
                cadat.DATA.Zl = sort(-cadat.DATA.Zl,'ascend');
                cadat.DATA.Zlc = sort(-cadat.DATA.Zlc,'ascend');
                cadat.DATA.Z = -cadat.DATA.Z;
                cadat.DATA.w = -cadat.DATA.w;
                flip_field_z = 1;
                flip_stl_z = 1;
        end
        % write STL x, y, z according to flipped STL data
        cadat.STL.x = stl_rm_x';
        cadat.STL.y = stl_rm_y';
        cadat.STL.zadj = stl_rm_z';
        cadat.STL.zplanes_adj = unique(cadat.STL.zadj);
        
        if cadat.OPTIONS.save_stl_flips
            % Save new stl file that was adjusted and flipped
            save([cadat.DIR.savefiles,cadat.STL.adjstlflipfile,'.mat'],'stl_rm_x','stl_rm_y','stl_rm_z')
        else
            % Load the STL flipping data
            fp = load([cadat.DIR.savefiles,cadat.STL.adjstlflipfile,'.mat']);
            cadat.STL.x = fp.stl_rm_x';
            cadat.STL.y = fp.stl_rm_y';
            cadat.STL.zadj = fp.stl_rm_z';
            cadat.STL.zplanes_adj = unique(cadat.STL.zadj);
        end
        
        % Save if data was flipped
        save([cadat.DIR.savefiles,cadat.DATA.flipfilename],'flip_field_x','flip_field_y','flip_field_z')
        fprintf('\n...completed successfully');
    else % If flipping has already been determined
        % Load adjusted and flipped STL file
        fprintf('\nImporting STL and flipping data ...')
        fp = load([cadat.DIR.savefiles,cadat.STL.adjstlflipfile,'.mat']);
        cadat.STL.x = fp.stl_rm_x';
        cadat.STL.y = fp.stl_rm_y';
        cadat.STL.zadj = fp.stl_rm_z';
        cadat.STL.zplanes_adj = unique(cadat.STL.zadj);
        
        % Load data flipping needed
        fp = load([cadat.DIR.savefiles,cadat.DATA.flipfilename]);
        cadat.DATA.flipfieldX = fp.flip_field_x;
        cadat.DATA.flipfieldY = fp.flip_field_y;
        cadat.DATA.flipfieldZ = fp.flip_field_z;
        % Flip x and u, if needed
        if fp.flip_field_x
            cadat.DATA.x = -cadat.DATA.x;
            cadat.DATA.xf = -cadat.DATA.xf;
            cadat.DATA.Xl = sort(-cadat.DATA.Xl,'ascend');
            cadat.DATA.Xlc = sort(-cadat.DATA.Xlc,'ascend');
            cadat.DATA.X = -cadat.DATA.X;
            cadat.DATA.u = -cadat.DATA.u;
            cadat.DATA.centvalX = -cadat.DATA.centvalX;
        end
        % Flip y and v, if needed
        if fp.flip_field_y
            cadat.DATA.y = -cadat.DATA.y;
            cadat.DATA.yf = -cadat.DATA.yf;
            cadat.DATA.Yl = sort(-cadat.DATA.Yl,'ascend');
            cadat.DATA.Ylc = sort(-cadat.DATA.Ylc,'ascend');
            cadat.DATA.Y = -cadat.DATA.Y;
            cadat.DATA.v = -cadat.DATA.v;
            cadat.DATA.centvalY = -cadat.DATA.centvalY;
        end
        % Flip z and w, if needed
        if fp.flip_field_z
            cadat.DATA.z = -cadat.DATA.z;
            cadat.DATA.zf = -cadat.DATA.zf;
            cadat.DATA.Zl = sort(-cadat.DATA.Zl,'ascend');
            cadat.DATA.Zlc = sort(-cadat.DATA.Zlc,'ascend');
            cadat.DATA.Z = -cadat.DATA.Z;
            cadat.DATA.w = -cadat.DATA.w;
            cadat.DATA.centvalZ = -cadat.DATA.centvalZ;
        end
        fprintf('...completed successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CREATE MASK BY CREATING CLOSED CONTOURS FROM STL DATA POINTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This step creates closed contour masks using the adjusted STL files.
    % To do this, the user needs to accept each mask, or choose to adjust
    % it by hand.
    if (cadat.OPTIONS.CREATEMASK) &&  (vnum == cadat.TIME.ts) && (cadat.OPTIONS.TESTRESOLUTION == 0)
        fprintf('\nCreating mask from STL...\n')
        % Initialize mask
        x_m = (min(cadat.STL.x)-cadat.STL.dx*8):cadat.STL.dx:(max(cadat.STL.x)+cadat.STL.dx*8);
        y_m = (min(cadat.STL.y)-cadat.STL.dy*8):cadat.STL.dy:(max(cadat.STL.y)+cadat.STL.dy*8);
        cadat.STL.x_m = x_m;
        cadat.STL.y_m = y_m;
        
        iter1prompt = 'Accept Mask? [1-Yes, 0-No]: ';
        ccprompt = 'Adjusted connected segments? [0-No, -1-Update Permenantly ,Else-new threshold]: ';
        ccprompt_perm = 'Adjust connected segments permenantly to [new threshold]: ';
        prompt = 'Accept Mask? [1-Yes, 0-No, -1-Move on, -2-Use Boundary]: ';
        adjconsol_prompt = 'Update consolidation threshold distance? [0-No, Else-new threshold]: ';
        adjneigh_prompt = 'Update neighbor threshold distance? [0-No, Else-new threshold]: ';
        promptadjcurrvals = 'Update consol/neigh thresholds permenantly? [1-Yes, 0-No]: ';
        consol_thresh_curr = 0.5;
        neigh_thresh_curr = 4;
        conn_pts_thresh = 2;
        for zz = 1:1:length(cadat.STL.zplanes_adj)
            % Get points in current z-plane
            stl_z_pts = find(cadat.STL.zadj == cadat.STL.zplanes_adj(zz));
            stl_plane = [cadat.STL.x(stl_z_pts),cadat.STL.y(stl_z_pts)]; % X and Y points of plane

            % Identify connected points within the chosen z-plane in order to create
            % closed contour boundaries - first iteration
            [pt_segs] = connected_points_v2(stl_plane(:,1),stl_plane(:,2),conn_pts_thresh);
            
            mask_cz = zeros(length(y_m),length(x_m)); % Initialize mask plane
            num_segs = length(unique(pt_segs));
            for sn = 1:1:num_segs
                % Initialize number of iterations
                iter = 1;
                % Determine wheter or not to use boundary
                num_pts_seg = sum(pt_segs == sn);
                use_bnd_method = num_pts_seg < 250;
                % Initialize thresholds
                consol_thresh = consol_thresh_curr;
                neigh_thresh = neigh_thresh_curr;
                % Get mask from ordered points approach
                if use_bnd_method == 0
                    [mask_init] = create_bnds_CCpoints_v3(stl_plane(:,1),stl_plane(:,2),pt_segs,sn,x_m,y_m,iter,consol_thresh,neigh_thresh);
                else
                    % Use a simple boundary in this case
                    x_pts = stl_plane(:,1);
                    y_pts = stl_plane(:,2);
                    % Only keep points in current segment
                    segi = find(pt_segs == sn);    
                    x_pts = x_pts(segi);
                    y_pts = y_pts(segi);

                    % Convert x,y from mm to pixels
                    ox_pix = zeros(size(x_pts));
                    oy_pix = zeros(size(y_pts));
                    for q = 1:1:length(x_pts)
                        % Get current x, y point
                        curr_x = x_pts(q);
                        curr_y = y_pts(q);

                        % Get closest index in mask
                        [~,x_ind] = min(abs(x_m-curr_x));
                        [~,y_ind] = min(abs(y_m-curr_y));

                        % Save closest index
                        ox_pix(q) = x_ind;
                        oy_pix(q) = y_ind;
                    end

                    % Compute boundary
                    kbnd = boundary(ox_pix,oy_pix);

                    % Convert boundary to mask
                    mask_init = poly2mask(ox_pix(kbnd),oy_pix(kbnd),length(y_m),length(x_m));
                end

                % Smooth the mask
                mask_init_sm = wiener2(mask_init,[5 5]);
                mask_init_sm = mask_init_sm > 0;

                % Add new mask to current z-plane mask
                mask_cz = mask_cz + mask_init_sm;

            end
            
            % Determine if first iteration is acceptable
            figure(15); hold off; imagesc(x_m,y_m,mask_cz); colormap gray
            set(gca,'FontSize',14)
            title(['Segment Mask for Z Index Plane = ',num2str(zz)])
            hold on; scatter(stl_plane(:,1),stl_plane(:,2),4,'r')
            user_reject = input(iter1prompt);
            
            %%% Re-iterate with more choices if the user rejects %%%%%%%%%%
            if user_reject ~= 1
                % Check if the connected points is right, if not have user
                % adjust to combine connected segments
                CC_okay = 0;
                num_segs = length(unique(pt_segs));
                x_all = stl_plane(:,1);
                y_all = stl_plane(:,2);
                figure(20); hold off; axis([-35 35 -35 35])
                for sn = 1:1:num_segs
                    segi = find(pt_segs == sn);
                    figure(20); plot(x_all(segi),y_all(segi),'o'); hold on;
                end
                set(gca,'FontSize',14); title(['Connected Segments for Z Index Plane = ',num2str(zz)])
                while CC_okay == 0
                    pt_dist_thresh = input(ccprompt);
                    if pt_dist_thresh == 0
                        CC_okay = 1;
                    elseif pt_dist_thresh == -1
                        conn_pts_thresh = input(ccprompt_perm);
                        pt_dist_thresh = conn_pts_thresh;
                        [pt_segs] = connected_points_v2(stl_plane(:,1),stl_plane(:,2),pt_dist_thresh);
                        num_segs = length(unique(pt_segs));
                        x_all = stl_plane(:,1);
                        y_all = stl_plane(:,2);
                        figure(20); hold off;
                        for sn = 1:1:num_segs
                            segi = find(pt_segs == sn);
                            figure(20); plot(x_all(segi),y_all(segi),'o'); hold on;
                        end
                        axis([-35 35 -35 35])
                        set(gca,'FontSize',14); title(['Connected Segments for Z Index Plane = ',num2str(zz)])
                    else
                        [pt_segs] = connected_points_v2(stl_plane(:,1),stl_plane(:,2),pt_dist_thresh);
                        num_segs = length(unique(pt_segs));
                        x_all = stl_plane(:,1);
                        y_all = stl_plane(:,2);
                        figure(20); hold off;
                        for sn = 1:1:num_segs
                            segi = find(pt_segs == sn);
                            figure(20); plot(x_all(segi),y_all(segi),'o'); hold on;
                        end
                        axis([-35 35 -35 35])
                        set(gca,'FontSize',14); title(['Connected Segments for Z Index Plane = ',num2str(zz)])
                    end
                end

                % Create closed boundary contours using the point segments identified
                mask_cz = zeros(length(y_m),length(x_m)); % Initialize mask plane

                % Iterate through all point segments to create mask. For first
                % iteration, algorithm will process all and if user accepts, it
                % moves on. Otherwise it goes through all steps
                num_segs = length(unique(pt_segs));
                for sn = 1:1:num_segs
                    % Initialize number of iterations to 2 since first pass
                    % was rejected
                    iter = 2;
                    % Initialize user prompt (user response) variable
                    user_prompt = 0;
                    % Initialize thresholds
                    consol_thresh = consol_thresh_curr;
                    neigh_thresh = neigh_thresh_curr;
                    while user_prompt ~= 1
                        % Get mask from ordered points approach
                        if user_prompt ~= -2
                            [mask_init] = create_bnds_CCpoints_v3(stl_plane(:,1),stl_plane(:,2),pt_segs,sn,x_m,y_m,iter,consol_thresh,neigh_thresh);
                        else
                            % Use a simple boundary in this case
                            x_pts = stl_plane(:,1);
                            y_pts = stl_plane(:,2);
                            % Only keep points in current segment
                            segi = find(pt_segs == sn);    
                            x_pts = x_pts(segi);
                            y_pts = y_pts(segi);

                            % Convert x,y from mm to pixels
                            ox_pix = zeros(size(x_pts));
                            oy_pix = zeros(size(y_pts));
                            for q = 1:1:length(x_pts)
                                % Get current x, y point
                                curr_x = x_pts(q);
                                curr_y = y_pts(q);

                                % Get closest index in mask
                                [~,x_ind] = min(abs(x_m-curr_x));
                                [~,y_ind] = min(abs(y_m-curr_y));

                                % Save closest index
                                ox_pix(q) = x_ind;
                                oy_pix(q) = y_ind;
                            end

                            % Compute boundary
                            kbnd = boundary(ox_pix,oy_pix);

                            % Convert boundary to mask
                            mask_init = poly2mask(ox_pix(kbnd),oy_pix(kbnd),length(y_m),length(x_m));
                        end

                        % Smooth the mask
                        mask_init_sm = wiener2(mask_init,[5 5]);
                        mask_init_sm = double(mask_init_sm > 0); % Ensure mask is only 1's and 0's

                        % Determine if mask is acceptable, if not re-run with different
                        % distance thresholds or starting points (all user
                        % selected)
                        figure(15); hold off; imagesc(x_m,y_m,mask_init_sm); colormap gray
                        set(gca,'FontSize',14)
                        title(['Segment Mask for Z Index Plane = ',num2str(zz)])
                        hold on; scatter(stl_plane(:,1),stl_plane(:,2),4,'r')
                        user_prompt = input(prompt);
                        if user_prompt == 0 % If user rejects, see if consolidation, neighbor thresholds need to be adjusted
                            update_permanent = input(promptadjcurrvals);
                            new_thresh = input(adjconsol_prompt);
                            if new_thresh ~= 0, consol_thresh = new_thresh; end
                            new_thresh = input(adjneigh_prompt);
                            if new_thresh ~= 0, neigh_thresh = new_thresh; end
                            if update_permanent
                                consol_thresh_curr = consol_thresh;
                                neigh_thresh_curr = neigh_thresh;
                            end
                        elseif user_prompt == -1
                            user_prompt = 1;
                            fprintf('User unhappy - move on\n!');
                        end

                        % Number of iterations on current mask portion
                        iter = iter + 1;
                    end

                    % Add new mask to current z-plane mask
                    mask_cz = mask_cz + mask_init_sm;
                end
            end
            % Ensure mask is between zero and 1
            mask_cz = mask_cz - min(mask_cz(:));
            mask_cz = mask_cz/max(mask_cz(:));
            
            % Save final mask
            mask(:,:,zz) = mask_cz;

            % Plot and print the final mask 
            figure(25); hold off; imagesc(x_m,y_m,mask(:,:,zz)); colormap gray;
            set(gca,'FontSize',14); set(gcf,'Color',[1 1 1]);
            title(['Final Mask for Z index ',num2str(zz)])
            hold on; scatter(stl_plane(:,1),stl_plane(:,2),4,'r')
            drawnow;
            %print(25,'-djpeg',['~/Desktop/ca_pics/stl_masks_v2/STL_pts_',num2str(zz,'%02i'),'.jpg'])
        end
        % Save final mask
        cadat.STL.mask = mask;
        save([cadat.DIR.savefiles,cadat.STL.maskfilename],'mask','x_m','y_m')
        fprintf('\n...mask created successfully\n');
    else
        fprintf('\nImporting STL mask...')
        fp = load([cadat.DIR.savefiles,cadat.STL.maskfilename]);
        cadat.STL.mask = fp.mask;
        cadat.STL.x_m = fp.x_m;
        cadat.STL.y_m = fp.y_m;
        fprintf('completed successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CREATE DATA MASK OR IMPORT IT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Create mask of points of non-zero velocity in field %%%%%%%%%%%%%%%
    if ((cadat.OPTIONS.CREATEDATAMASK) &&  (vnum == cadat.TIME.ts)) || (cadat.OPTIONS.TESTRESOLUTION == 1)
        if cadat.OPTIONS.TESTRESOLUTION == 0
            fprintf('\nCreating data velocity mask...')
        else
            fprintf('\nCreating data velocity mask [RESOLUTION TEST ONLY]...')
        end
        % Initialize velocity mask
        cadat.DATA.VelMask = zeros(size(cadat.DATA.U));
        if strcmp(cadat.FILE.filetype,'TOMO')
            % For tomo, only need one file to load in (already loaded)
            for zz = 1:1:size(cadat.DATA.U,3)
                cadat.DATA.VelMask(:,:,zz) = abs(cadat.DATA.U(:,:,zz)) > 0;
            end
        else % For STB, CFD, and MRI, want to use all velocity files
            prev_mark = 0;
            itstep_ct = 1;
            for itstep = cadat.VELMASK.ts:cadat.VELMASK.tstep:cadat.VELMASK.te
                %%% Import velocity file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(cadat.FILE.filetype,'STB')
                    fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(itstep,'%06i'),'.mat']);
                else
                    fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(itstep,'%05i'),'.mat']);
                end
                % Original x,y,z coordinates need to be adjusted as they are in
                % original import step
                % For MRI, need to convert to mm again and remove zero
                % points again
                if cadat.DATA.inMM == 0
                    cx = fp.x * 1000; % Convert to mm
                    cy = fp.y * 1000;
                    cz = fp.z * 1000;
                else
                    cx = fp.x;
                    cy = fp.y;
                    cz = fp.z;
                end
                % Keep only the data that has non-zero velocities
                vel_kp = (abs(fp.u) + abs(fp.v) + abs(fp.w)) > 0;
                cx = cx(vel_kp == 1);
                cy = cy(vel_kp == 1);
                cz = cz(vel_kp == 1);
                cu = fp.u(vel_kp == 1);
                cv = fp.v(vel_kp == 1);
                cw = fp.w(vel_kp == 1); 
                
                % Need to flip the data if data flipping was done. Only
                % need to flip the coordinates, the magnitude of the
                % velocities is not important
                % Flip x, if needed
                if cadat.DATA.flipfieldX
                    cx = -cx;
                end
                % Flip y, if needed
                if cadat.DATA.flipfieldY
                    cy = -cy;
                end
                % Flip z and w, if needed
                if cadat.DATA.flipfieldZ
                    cz = -cz;
                end  
                
                % Need to center
                cx = cx - cadat.DATA.centvalX;
                cy = cy - cadat.DATA.centvalY;
                cz = cz - cadat.DATA.centvalZ;

                % For each velocity point, add it to the mask location
                for zz = 1:1:size(cx)
                    % Get current x,y,z points
                    curr_x = cx(zz);
                    curr_y = cy(zz);
                    curr_z = cz(zz);

                    % Find the closest x, y, z coordinate
                    [~,r_ind] = min(abs(cadat.DATA.Ylc - curr_y));
                    [~,c_ind] = min(abs(cadat.DATA.Xlc - curr_x));
                    [~,z_ind] = min(abs(cadat.DATA.Zlc - curr_z));
                    cadat.DATA.VelMask(r_ind,c_ind,z_ind) = 1;
                end
                % Iterate the count
                itstep_ct = itstep_ct + 1;
                
                % Print the percent completed for steps of 10
                perc_complete = floor(100*itstep_ct/length(cadat.VELMASK.ts:cadat.VELMASK.tstep:cadat.VELMASK.te));
                if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                    fprintf(' %i%%',perc_complete)
                    prev_mark = perc_complete;
                end
                
            end
        end
        
        % If resolution is being tested, program will terminate after
        % plotting the data mask
        if cadat.OPTIONS.TESTRESOLUTION == 1
            mid_z = round(length(cadat.DATA.Zlc)/2);
            close all; figure(20); imagesc(cadat.DATA.VelMask(:,:,mid_z));
            set(gca,'FontSize',14); title(['Data Mask Z-Plane = ',num2str(mid_z)])
            return
        else
            %%% Adjust resolution to be same as STL mask %%%%%%%%%%%%%%%%%%
            % If the resolution of the data is not the same as that of the
            % STL, the data mask needs to be super-sampled
            if cadat.DATA.dx ~= cadat.STL.dx
                spmultfact = round(cadat.DATA.dx/cadat.STL.dx);
                
                % Need to add enough rows, columns, z-planes to make the
                % velocity mask the same resolution as the STL mask
                totR = size(cadat.DATA.VelMask,1)*spmultfact;
                totC = size(cadat.DATA.VelMask,2)*spmultfact;
                totZ = size(cadat.DATA.VelMask,3)*spmultfact;
                
                % Save the original mask
                cadat.DATA.VelMask_orig = cadat.DATA.VelMask;
                
                % Create new mask and coordinates
                new_x = min(cadat.DATA.Xlc):cadat.STL.dx:max(cadat.DATA.Xlc);
                new_y = min(cadat.DATA.Ylc):cadat.STL.dy:max(cadat.DATA.Ylc);
                new_z = min(cadat.DATA.Zlc):cadat.STL.dz:max(cadat.DATA.Zlc);
                new_mask = zeros(length(new_y),length(new_x),length(new_z));
                mask_try = imresize(cadat.DATA.VelMask,[length(new_y),length(new_x)]);
                mask_try = double(mask_try > 0.03);
                mask_try2 = permute(mask_try,[2 3 1]);
                mask_try2 = imresize(mask_try2,[length(new_x),length(new_z)]);
                mask_try2 = permute(mask_try2,[3 1 2]);
                mask_tryf = double(mask_try2 > 0.03);
                cadat.DATA.VelMask = mask_tryf;
%                 %%% Add the original mask to the new mask %%%%%%%%%%%%%%%%%
%                 % Create lines of the mask and X Y Z coordinates
%                 velmaskI = reshape(cadat.DATA.VelMask,[size(cadat.DATA.VelMask,1)*size(cadat.DATA.VelMask,2)*size(cadat.DATA.VelMask,3),1]);
%                 [mX,mY,mZ] = meshgrid(cadat.DATA.Xlc,cadat.DATA.Ylc,cadat.DATA.Zlc);
%                 velmaskX = reshape(mX,[size(cadat.DATA.VelMask,1)*size(cadat.DATA.VelMask,2)*size(cadat.DATA.VelMask,3),1]);
%                 velmaskY = reshape(mY,[size(cadat.DATA.VelMask,1)*size(cadat.DATA.VelMask,2)*size(cadat.DATA.VelMask,3),1]);
%                 velmaskZ = reshape(mZ,[size(cadat.DATA.VelMask,1)*size(cadat.DATA.VelMask,2)*size(cadat.DATA.VelMask,3),1]);
%                 
%                 % Keep only the non-zero points
%                 kp_pts = velmaskI > 0;
%                 velmaskI = velmaskI(kp_pts == 1);
%                 velmaskX = velmaskX(kp_pts == 1);
%                 velmaskY = velmaskY(kp_pts == 1);
%                 velmaskZ = velmaskZ(kp_pts == 1);
%                 
%                 % Iterate through all non-zero points and put them into the
%                 % resolution adjusted mask
%                 for ptn = 1:1:length(velmaskI)
%                     % Get current point
%                     curr_x = velmaskX(ptn);
%                     curr_y = velmaskY(ptn);
%                     curr_z = velmaskZ(ptn);
%                     
%                     % Find the point on the new mask
%                     rind = find(curr_y == new_y);
%                     cind = find(curr_x == new_x);
%                     zind = find(curr_z == new_z);
%                     
%                     % Place the point on the new index
%                     new_mask(rind,cind,zind) = 1;
%                 end
            end
            
            %%% Adjust data mask to be of same size as STL mask %%%%%%%%%%%
            % Adjust number of Z planes
            add_z = size(cadat.DATA.VelMask,3) - size(cadat.STL.mask,3); % Number of planes to add to mask
            add_plane = zeros(size(cadat.STL.mask(:,:,1)));
            add_planeD = zeros(size(cadat.DATA.VelMask(:,:,1)));
            if add_z > 0 % Need to add planes to STL mask, need to adjust the STL mask coordinates
                fnumplane = floor(add_z/2);
                bnumplane = ceil(add_z/2);
                fadd = repmat(add_plane,[1 1 fnumplane]);
                badd = repmat(add_plane,[1 1 bnumplane]);
                cadat.STL.maskz = cat(3,fadd,cadat.STL.mask,badd);
                cadat.DATA.VelMaskz = cadat.DATA.VelMask;
                fz_add = -cadat.STL.dz*(fnumplane:-1:1) + cadat.STL.zplanes_adj(1);
                bz_add = cadat.STL.dz*(1:1:bnumplane) + cadat.STL.zplanes_adj(end);
                stl_z = [fz_add';cadat.STL.zplanes_adj;bz_add'];
            elseif add_z < 0 % Need to add planes to velocity field mask
                fnumplane = floor(abs(add_z)/2);
                bnumplane = ceil(abs(add_z)/2);
                fadd = repmat(add_planeD,[1 1 fnumplane]);
                badd = repmat(add_planeD,[1 1 bnumplane]);
                cadat.DATA.VelMaskz = cat(3,fadd,cadat.DATA.VelMask,badd);
                cadat.STL.maskz = cadat.STL.mask;
                stl_z = cadat.STL.zplanes_adj;
            else
                cadat.STL.maskz = cadat.STL.mask;
                cadat.DATA.VelMaskz = cadat.DATA.VelMask;
                stl_z = cadat.STL.zplanes_adj;
            end
            % Adjust number of rows
            add_r = size(cadat.DATA.VelMask,1) - size(cadat.STL.mask,1);
            add_rows = zeros(1,size(cadat.STL.maskz,2),size(cadat.STL.maskz,3));
            add_rowsD = zeros(1,size(cadat.DATA.VelMaskz,2),size(cadat.DATA.VelMaskz,3));
            if add_r > 0 % Need to add more rows to STL mask, need to adjust the STL mask coordinates
                fnumrows = floor(abs(add_r)/2);
                bnumrows = ceil(abs(add_r)/2);
                fadd = repmat(add_rows,[fnumrows 1 1]);
                badd = repmat(add_rows,[bnumrows 1 1]);
                cadat.STL.maskyz = cat(1,fadd,cadat.STL.maskz,badd);
                cadat.DATA.VelMaskyz = cadat.DATA.VelMaskz;
                fy_add = -cadat.STL.dy*(fnumrows:-1:1) + cadat.STL.y_m(1);
                by_add = cadat.STL.dy*(1:1:bnumrows) + cadat.STL.y_m(end);
                stl_y = [fy_add,cadat.STL.y_m,by_add];
            elseif add_r < 0 % Need to add more rows to velocity field mask
                fnumrows = floor(abs(add_r)/2);
                bnumrows = ceil(abs(add_r)/2);
                fadd = repmat(add_rowsD,[fnumrows 1 1]);
                badd = repmat(add_rowsD,[bnumrows 1 1]);
                cadat.DATA.VelMaskyz = cat(1,fadd,cadat.DATA.VelMaskz,badd);
                cadat.STL.maskyz = cadat.STL.maskz;
                stl_y = cadat.STL.y_m;
            else
                cadat.STL.maskyz = cadat.STL.maskz;
                cadat.DATA.VelMaskyz = cadat.DATA.VelMaskz;
                stl_y = cadat.STL.y_m;
            end
            % Adjust number of columns
            add_c = size(cadat.DATA.VelMask,2) - size(cadat.STL.mask,2);
            add_cols = zeros(size(cadat.STL.maskyz(:,1,:)));
            add_colsD = zeros(size(cadat.DATA.VelMaskyz(:,1,:)));
            if add_c > 0 % Need to add more cols to STL mask
                fnumrows = floor(abs(add_c)/2);
                bnumrows = ceil(abs(add_c)/2);
                fadd = repmat(add_cols,[1 fnumrows 1]);
                badd = repmat(add_cols,[1 bnumrows 1]);
                cadat.STL.maskf = cat(2,fadd,cadat.STL.maskyz,badd);
                cadat.DATA.VelMaskf = cadat.DATA.VelMaskyz;
                fx_add = -cadat.STL.dx*(fnumrows:-1:1) + cadat.STL.x_m(1);
                bx_add = cadat.STL.dx*(1:1:bnumrows) + cadat.STL.x_m(end);
                stl_x = [fx_add,cadat.STL.x_m,bx_add];
            elseif add_c < 0 % Need to add more cols to velocity field mask
                fnumrows = floor(abs(add_c)/2);
                bnumrows = ceil(abs(add_c)/2);
                fadd = repmat(add_colsD,[1 fnumrows 1]);
                badd = repmat(add_colsD,[1 bnumrows 1]);
                cadat.DATA.VelMaskf = cat(2,fadd,cadat.DATA.VelMaskyz,badd);
                cadat.STL.maskf = cadat.STL.maskyz;
                stl_x = cadat.STL.x_m;
            else
                cadat.STL.maskf = cadat.STL.maskyz;
                cadat.DATA.VelMaskf = cadat.DATA.VelMaskyz;
                stl_x = cadat.STL.x_m;
            end

            % Save filtered masks
            VelMaskf = cadat.DATA.VelMaskf;
            stlmaskf = cadat.STL.maskf;
            save([cadat.DIR.savefiles,cadat.VELMASK.filtmaskname],'VelMaskf')
            save([cadat.DIR.savefiles,cadat.STL.filtmaskname],'stlmaskf','stl_x','stl_y','stl_z')
        end
    else
        fprintf('\nImporting Data mask...')
        % Load filtered masks
        fp = load([cadat.DIR.savefiles,cadat.VELMASK.filtmaskname]);
        cadat.DATA.VelMaskf = fp.VelMaskf;
        fp = load([cadat.DIR.savefiles,cadat.STL.filtmaskname]);
        cadat.STL.maskf = fp.stlmaskf;
        stl_x = fp.stl_x;
        stl_y = fp.stl_y;
        stl_z = fp.stl_z;
        fprintf('completed successfully')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIND OR IMPORT OPTIMAL SHIFT FOR REGISTRATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (cadat.OPTIONS.FINDOPTSHIFT) &&  (vnum == cadat.TIME.ts)
        fprintf('\nFinding optimal registration shift...\n')
        %%% Find x, y, z displacements that maximize point overlap between
        %%% the STL and velocity masks. Start by having user selected a
        %%% guess for what the offsets should be to minimize searching
        %%% needed.
        % Plot the middle z-index for X,Y user guesses
        plt_z_plane = round(size(cadat.STL.maskf,3)/2);
        figure(60); subplot(1,2,1); imagesc(cadat.STL.maskf(:,:,plt_z_plane));
        colormap gray; set(gca,'FontSize',14); title('STL Mask')
        subplot(1,2,2); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane));
        colormap gray; set(gca,'FontSize',14); title('Data Mask'); set(gcf,'Color',[1 1 1]);
        prompt_x = 'X-Offset Guess [Side-Side, + Right, - Left]:  ';
        prompt_y = 'Y-Offset Guess [Up-Down, + Down, - Up]:  ';
        prompt_z = 'Z-Offset Guess [Forward-Back, Pick title]:  ';
        x_guess = input(prompt_x);
        y_guess = input(prompt_y);
        % Plot multiple z planes for Z user guess
        figure(60); subplot(3,3,1); imagesc(cadat.STL.maskf(:,:,plt_z_plane)); colormap gray;
        set(gca,'FontSize',14); title('STL -  Z = 0')
        subplot(3,3,2); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = 0')
        subplot(3,3,3); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane+2)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = -2')
        subplot(3,3,4); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane+4)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = -4')
        subplot(3,3,5); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane+6)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = -6')
        subplot(3,3,6); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane+8)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = -8')
        subplot(3,3,7); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane-2)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = +2')
        subplot(3,3,8); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane-4)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = +4')
        subplot(3,3,9); imagesc(cadat.DATA.VelMaskf(:,:,plt_z_plane-6)); colormap gray; set(gca,'FontSize',14); title('Data -  Z = +6')
        z_guess = input(prompt_z);
        
        % Set search vectors based on guesses
        x_srch = x_guess-5:1:x_guess+5;
        y_srch = y_guess-5:1:y_guess+5;
        z_srch = z_guess-5:1:z_guess+5;
        
        % Initialize output and counters
        srch_pts = zeros(length(y_srch),length(x_srch),length(z_srch));
        num_srch_pts = size(srch_pts,1)*size(srch_pts,2)*size(srch_pts,3);
        iter = 1;
        prev_mark = 0;
        %%% FIND OPTIMAL TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Searching using guess (%i, %i, %i)...',x_guess, y_guess, z_guess)
        tic
        for xd = 1:1:length(x_srch)
            for yd = 1:1:length(y_srch)
                for zd = 1:1:length(z_srch)
                    % Get current index shifts
                    xshift = x_srch(xd);
                    yshift = y_srch(yd);
                    zshift = z_srch(zd);

                    % Place the velocity fields on the zero fields based on the shifts
                    % Complete z-shift
                    if zshift > 0
                        zadd = zeros(size(cadat.DATA.VelMaskf,1),size(cadat.DATA.VelMaskf,2),zshift);
                        vel_zshift = cat(3,zadd,cadat.DATA.VelMaskf(:,:,1:end-zshift));
                    elseif zshift < 0
                        zadd = zeros(size(cadat.DATA.VelMaskf,1),size(cadat.DATA.VelMaskf,2),abs(zshift));
                        vel_zshift = cat(3,cadat.DATA.VelMaskf(:,:,abs(zshift)+1:end),zadd);
                    else
                        vel_zshift = cadat.DATA.VelMaskf;
                    end

                    % Complete x-shift
                    % Adjust columns based on x shift value
                    if xshift > 0
                        xadd = zeros(size(vel_zshift,1),xshift,size(vel_zshift,3));
                        vel_xshift = cat(2,xadd,vel_zshift(:,1:end-xshift,:));
                    elseif xshift < 0
                        xadd = zeros(size(vel_zshift,1),abs(xshift),size(vel_zshift,3));
                        vel_xshift = cat(2,vel_zshift(:,abs(xshift-1):end,:),xadd);
                    else
                        vel_xshift = vel_zshift;
                    end

                    % Complete y-shift
                    % Adjust rows based on y shift value
                    if yshift > 0
                        yadd = zeros(yshift,size(vel_xshift,2),size(vel_xshift,3));
                        vel_shifted = cat(1,yadd,vel_xshift(1:end-yshift,:,:));
                    elseif yshift < 0
                        yadd = zeros(abs(yshift),size(vel_xshift,2),size(vel_xshift,3));
                        vel_shifted = cat(1,vel_xshift(abs(yshift-1):end,:,:),yadd);
                    else
                        vel_shifted = vel_xshift;
                    end
                    
                    % Determine number of overlapping points at displacement
                    num_non0_pts = vel_shifted.*cadat.STL.maskf;
                    srch_pts(yd,xd,zd) = sum(num_non0_pts(:))/(size(vel_shifted,1)*size(vel_shifted,2)*size(vel_shifted,3));

                    % Print the percent completed for steps of 10
                    perc_complete = floor(100*iter/num_srch_pts);
                    if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                        fprintf(' %i%%',perc_complete)
                        prev_mark = perc_complete;
                    end
                    tval = toc;
                    if iter == 1
                        testimate = round(num_srch_pts*tval);
                        testimate_hr = floor(testimate/(60*60));
                        testimate_min = floor((testimate-testimate_hr*60*60)/60);
                        testimate_sec = testimate - (testimate_hr*60*60 + testimate_min*60);
                        fprintf('\nEstimated time to complete: %02i:%02i:%02i \n',testimate_hr,testimate_min,testimate_sec)
                        fprintf('Percent Completed: ')
                    end
                    iter = iter + 1;
                end
            end
        end
        
        % Locate optimal shift (maximum overlap)
        [msrch,mzind] = max(srch_pts,[],3);
        [msrch,mxind] = max(msrch,[],2);
        [msrch,myind] = max(msrch);
        y_opt = myind;
        x_opt = mxind(myind);
        z_opt = mzind(y_opt,x_opt);
        % Collect the center values
        centvalX = cadat.DATA.centvalX;
        centvalY = cadat.DATA.centvalY;
        centvalZ = cadat.DATA.centvalZ;
        
        % Save optimal shift as well as the center values used for it
        save([cadat.DIR.savefiles,cadat.SHIFT.optshiftfile],'x_opt','y_opt','z_opt','x_srch','y_srch','z_srch','centvalX','centvalY','centvalZ');
        
        % Save optimal shift in cadat structure for further calculations
        cadat.DATA.YRshift = y_srch(y_opt);
        cadat.DATA.XRshift = x_srch(x_opt);
        cadat.DATA.ZRshift = z_srch(z_opt);
        
        
        %%% FIND OPTIMAL ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rotation guesses
        rot_x_theta = 8;
        rot_y_theta = 0;
        rot_z_theta = 1;
        
        % Create X,Y,Z matrices
        [Xsr,Ysr,Zsr] = meshgrid(stl_x,stl_y,stl_z);
        % Reshape matrices into line
        Xsr_l = transpose(reshape(Xsr,[size(Xsr,1)*size(Xsr,2)*size(Xsr,3),1,1]));
        Ysr_l = transpose(reshape(Ysr,[size(Xsr,1)*size(Xsr,2)*size(Xsr,3),1,1]));
        Zsr_l = transpose(reshape(Zsr,[size(Xsr,1)*size(Xsr,2)*size(Xsr,3),1,1]));
        
        % Convert the x, y, z coordinates based on the shift,
        % centering, and conversion to stl coordinates
        % Center coordinates
        xc = Xsr_l;% - centvalX;
        yc = Ysr_l;% - centvalY;
        zc = Zsr_l;% - centvalZ;
        
        % Adjust based on shift
        Xsr_lc = xc + cadat.DATA.XRshift*cadat.STL.dx;
        Ysr_lc = yc + cadat.DATA.YRshift*cadat.STL.dy;
        Zsr_lc = zc + cadat.DATA.ZRshift*cadat.STL.dz;
        
        
        % Adjust the points based on the optimal mask
        % Reshape the mask into a line
        dmask_l = transpose(reshape(cadat.DATA.VelMaskf,[size(Xsr,1)*size(Xsr,2)*size(Xsr,3),1]));
        % Keep only the non-zero points
        kp_inds = (dmask_l>0);
        Xsr_lf = Xsr_l(kp_inds==1);
        Ysr_lf = Ysr_l(kp_inds==1);
        Zsr_lf = Zsr_l(kp_inds==1);
        Xsr_lcf = Xsr_lc(kp_inds==1);
        Ysr_lcf = Ysr_lc(kp_inds==1);
        Zsr_lcf = Zsr_lc(kp_inds==1);
        
        % Subsample the points to be around 2000
        sub_size = floor(length(Zsr_lf)/2000);
        Xsr_lfss = Xsr_lcf(1:sub_size:end);
        Ysr_lfss = Ysr_lcf(1:sub_size:end);
        Zsr_lfss = Zsr_lcf(1:sub_size:end);
        
        % Plot the STL and scatter plots
        stl_mask = cadat.STL.maskf;
        figure(32); close(gcf);
        figure(32);
        vv = double(stl_mask > 0);
        p2 = patch(isosurface(Xsr,Ysr,Zsr,vv,0));
        p2.FaceColor = [0.7 0.7 0.7];
        p2.EdgeColor = 'none';
        p2.FaceAlpha = 0.2;
        hold on; scatter3(Xsr_lf,Ysr_lf,Zsr_lf,3,'filled')
        scatter3(Xsr_lcf,Ysr_lcf,Zsr_lcf,3,'filled')
        pause(1E-3);
        
        % Need to center the XYZ points
        medX = median(Xsr_lfss);
        medY = median(Ysr_lfss);
        medZ = median(Zsr_lfss);
        Xsr_lfssc = Xsr_lfss - medX;
        Ysr_lfssc = Ysr_lfss - medY;
        Zsr_lfssc = Zsr_lfss - medZ;
        
        cXYZ = [Xsr_lfssc;Ysr_lfssc;Zsr_lfssc];
        % Initialze teh output
        rotSrch = zeros(length(rot_y_theta),length(rot_x_theta),length(rot_z_theta));
        iter = 1; prev_mark = 0;
        num_srch_pts = size(rotSrch,1)*size(rotSrch,2)*size(rotSrch,3);
        fprintf('\nSearching optimal rotation...')
        tic
        for qx = 1:1:length(rot_x_theta)
            for qy = 1:1:length(rot_y_theta)
                for qz = 1:1:length(rot_z_theta)
                    %%% Complete rotations of data %%%%%%%%%%%%%%%%%%%%%%%%
                    % X rotation
                    thetaX = rot_x_theta(qx);
                    % Get the X rotation matrix
                    rotmat = rotx(thetaX);
                    % Administer the X rotation
                    cXYZ_rotX = rotmat*cXYZ;

                    thetaY = rot_y_theta(qy);
                    % Get the X rotation matrix
                    rotmat = roty(thetaY);
                    % Administer the Y rotation
                    cXYZ_rotXY = rotmat*cXYZ_rotX;

                    thetaZ = rot_z_theta(qz);
                    % Get the Z rotation matrix
                    rotmat = rotz(thetaZ);
                    % Administer the Z rotation
                    cXYZ_rotXYZ = rotmat*cXYZ_rotXY;

                    % Find how many points are in the mask after the
                    % rotation
                    num_non_zero_pts = 0;
                    for cpt = 1:1:size(cXYZ_rotXYZ,2)
                        % Get the current X,Y,Z
                        currX = cXYZ_rotXYZ(1,cpt);
                        currY = cXYZ_rotXYZ(2,cpt);
                        currZ = cXYZ_rotXYZ(3,cpt);
                        % Find the nearest point on the
                        % rotated X,Y,Z coordinates
                        pt_dist = sqrt((Xsr-currX).^2 + (Ysr-currY).^2 + (Zsr-currZ).^2);
                        [minDist,minZ] = min(pt_dist,[],3);
                        [minMat,minX] = min(minDist,[],2);
                        [minPt,minY] = min(minMat);
                        % Get the minimum indices
                        indY = minY;
                        indX = minX(indY);
                        indZ = minZ(indY,indX);
                        % Place the correct value in
                        % the current index of the
                        % shifted matrix
                        inMask = cadat.STL.maskf(indY,indX,indZ);
                        if inMask
                            num_non_zero_pts = num_non_zero_pts + 1;
                        end
                    end
                    % Get a percent of non-zero points in mask
                    perc_non_zero = 100*(num_non_zero_pts/size(cXYZ_rotXYZ,2));
                    % Save the number of non-zero points
                    rotSrch(qy,qx,qz) = perc_non_zero;
                    
                    % Print the percent completed for steps of 10
                    perc_complete = floor(100*iter/num_srch_pts);
                    if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                        fprintf(' %i%%',perc_complete)
                        prev_mark = perc_complete;
                    end
                    if iter == 1
                        tval = toc;
                        testimate = round(num_srch_pts*tval);
                        testimate_hr = floor(testimate/(60*60));
                        testimate_min = floor((testimate-testimate_hr*60*60)/60);
                        testimate_sec = testimate - (testimate_hr*60*60 + testimate_min*60);
                        fprintf('\nEstimated time to complete: %02i:%02i:%02i \n',testimate_hr,testimate_min,testimate_sec)
                        fprintf('Percent Completed: ')
                    end
                    iter = iter + 1;
                    
                end
            end
        end
        [msrch,mqzind] = max(rotSrch,[],3);
        [msrch,mqxind] = max(msrch,[],2);
        [msrch,mqyind] = max(msrch);
        qy_opt = mqyind;
        qx_opt = mqxind(qy_opt);
        qz_opt = mqzind(qy_opt,qx_opt);
        
        rot_optX = 0;%rot_x_theta(qx_opt);
        rot_optY = 0;%rot_y_theta(qy_opt);
        rot_optZ = 0;%rot_z_theta(qz_opt);     
        
        % Plot the STL and scatter plots to ensure it is the correct
        % rotation
        % X rotation
        % Get the X rotation matrix
        rotmat = rotx(rot_optX);
        % Administer the X rotation
        cXYZ_rotX = rotmat*cXYZ;
        % Get the X rotation matrix
        rotmat = roty(rot_optY);
        % Administer the Y rotation
        cXYZ_rotXY = rotmat*cXYZ_rotX;
        % Get the Z rotation matrix
        rotmat = rotz(rot_optZ);
        % Administer the Z rotation
        cXYZ_rotXYZ = rotmat*cXYZ_rotXY;
        
        stl_mask = cadat.STL.maskf;
        figure(33); close(gcf);
        figure(33);
        vv = double(stl_mask > 0);
        p2 = patch(isosurface(Xsr,Ysr,Zsr,vv,0));
        p2.FaceColor = [0.7 0.7 0.7];
        p2.EdgeColor = 'none';
        p2.FaceAlpha = 0.2;
        hold on; %scatter3(Xsr_lfss,Ysr_lfss,Zsr_lfss,3,'filled');
        scatter3(cXYZ_rotXYZ(1,:)+medX,cXYZ_rotXYZ(2,:)+medY,cXYZ_rotXYZ(3,:)+medZ,3,'filled');
        pause(1E-3);

        save([cadat.DIR.savefiles,cadat.SHIFT.optrotfile],'rot_optX','rot_optY','rot_optZ','medX','medY','medZ');
        
        % Get the median and rotation values
        cadat.DATA.rot_optX = rot_optX;
        cadat.DATA.rot_optY = rot_optY;
        cadat.DATA.rot_optZ = rot_optZ;
        cadat.DATA.medX = medX;
        cadat.DATA.medY = medY;
        cadat.DATA.medZ = medZ;
    else
        fprintf('\nImporting optimal shift values...')
        fp = load([cadat.DIR.savefiles,cadat.SHIFT.optshiftfile]);
        cadat.DATA.XRshift = fp.x_srch(fp.x_opt);
        cadat.DATA.YRshift = fp.y_srch(fp.y_opt);
        cadat.DATA.ZRshift = fp.z_srch(fp.z_opt);
        % Get the center values
        cadat.DATA.centvalX = fp.centvalX;
        cadat.DATA.centvalY = fp.centvalY;
        cadat.DATA.centvalZ = fp.centvalZ;
        % Get the median and rotation values
        fp = load([cadat.DIR.savefiles,cadat.SHIFT.optrotfile]);
        cadat.DATA.rot_optX = fp.rot_optX;
        cadat.DATA.rot_optY = fp.rot_optY;
        cadat.DATA.rot_optZ = fp.rot_optZ;
        cadat.DATA.medX = fp.medX;
        cadat.DATA.medY = fp.medY;
        cadat.DATA.medZ = fp.medZ;
        fprintf('completed successfully')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REGISTER VELOCITY FIELD THEN MASK FROM ORIGINAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data velocity fields have already been flipped so only the
    % registration shifts need to be performed. Then the velocity field can
    % be masked.
    if cadat.OPTIONS.MASKREGVELOCITY
        fprintf('\nMasking and registering data velocity...') 
        %%% MASK DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % How to carry out this step changes based on the data being used.
        % TOMO uses a simple multiplication while STB, CFD, and MRI require
        % each point to be evaluated because this step is not intended to
        % also grid the data
        % Shift velocity based on optimal registration displacement
        % Get current index shifts
        xshift = cadat.DATA.XRshift;
        yshift = cadat.DATA.YRshift;
        zshift = cadat.DATA.ZRshift;
        if strcmp(cadat.FILE.filetype,'TOMO')
            %%% U SHIFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Place the velocity fields on the zero fields based on the shifts
            % Complete z-shift
            if zshift > 0
                zadd = zeros(size(cadat.DATA.U,1),size(cadat.DATA.U,2),zshift);
                vel_zshift = cat(3,zadd,cadat.DATA.U(:,:,1:end-zshift));
            elseif zshift < 0
                zadd = zeros(size(cadat.DATA.U,1),size(cadat.DATA.U,2),abs(zshift));
                vel_zshift = cat(3,cadat.DATA.U(:,:,abs(zshift)+1:end),zadd);
            else
                vel_zshift = cadat.DATA.U;
            end
            % Complete x-shift
            % Adjust columns based on x shift value
            if xshift > 0
                xadd = zeros(size(vel_zshift,1),xshift,size(vel_zshift,3));
                vel_xshift = cat(2,xadd,vel_zshift(:,1:end-xshift,:));
            elseif xshift < 0
                xadd = zeros(size(vel_zshift,1),abs(xshift),size(vel_zshift,3));
                vel_xshift = cat(2,vel_zshift(:,abs(xshift-1):end,:),xadd);
            else
                vel_xshift = vel_zshift;
            end
            % Complete y-shift
            % Adjust rows based on y shift value
            if yshift > 0
                yadd = zeros(yshift,size(vel_xshift,2),size(vel_xshift,3));
                vel_shifted = cat(1,yadd,vel_xshift(1:end-yshift,:,:));
            elseif yshift < 0
                yadd = zeros(abs(yshift),size(vel_xshift,2),size(vel_xshift,3));
                vel_shifted = cat(1,vel_xshift(abs(yshift-1):end,:,:),yadd);
            else
                vel_shifted = vel_xshift;
            end


            %%% V SHIFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Place the velocity fields on the zero fields based on the shifts
            % Complete z-shift
            if zshift > 0
                zadd = zeros(size(cadat.DATA.V,1),size(cadat.DATA.V,2),zshift);
                vel_zshift = cat(3,zadd,cadat.DATA.V(:,:,1:end-zshift));
            elseif zshift < 0
                zadd = zeros(size(cadat.DATA.V,1),size(cadat.DATA.V,2),abs(zshift));
                vel_zshift = cat(3,cadat.DATA.V(:,:,abs(zshift)+1:end),zadd);
            else
                vel_zshift = cadat.DATA.V;
            end
            % Complete x-shift
            % Adjust columns based on x shift value
            if xshift > 0
                xadd = zeros(size(vel_zshift,1),xshift,size(vel_zshift,3));
                vel_xshift = cat(2,xadd,vel_zshift(:,1:end-xshift,:));
            elseif xshift < 0
                xadd = zeros(size(vel_zshift,1),abs(xshift),size(vel_zshift,3));
                vel_xshift = cat(2,vel_zshift(:,abs(xshift-1):end,:),xadd);
            else
                vel_xshift = vel_zshift;
            end
            % Complete y-shift
            % Adjust rows based on y shift value
            if yshift > 0
                yadd = zeros(yshift,size(vel_xshift,2),size(vel_xshift,3));
                vel_shifted = cat(1,yadd,vel_xshift(1:end-yshift,:,:));
            elseif yshift < 0
                yadd = zeros(abs(yshift),size(vel_xshift,2),size(vel_xshift,3));
                vel_shifted = cat(1,vel_xshift(abs(yshift-1):end,:,:),yadd);
            else
                vel_shifted = vel_xshift;
            end
            cadat.DATA.Vr = vel_shifted.*cadat.STL.maskf;

            %%% W SHIFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Place the velocity fields on the zero fields based on the shifts
            % Complete z-shift
            if zshift > 0
                zadd = zeros(size(cadat.DATA.W,1),size(cadat.DATA.W,2),zshift);
                vel_zshift = cat(3,zadd,cadat.DATA.W(:,:,1:end-zshift));
            elseif zshift < 0
                zadd = zeros(size(cadat.DATA.W,1),size(cadat.DATA.W,2),abs(zshift));
                vel_zshift = cat(3,cadat.DATA.W(:,:,abs(zshift)+1:end),zadd);
            else
                vel_zshift = cadat.DATA.W;
            end
            % Complete x-shift
            % Adjust columns based on x shift value
            if xshift > 0
                xadd = zeros(size(vel_zshift,1),xshift,size(vel_zshift,3));
                vel_xshift = cat(2,xadd,vel_zshift(:,1:end-xshift,:));
            elseif xshift < 0
                xadd = zeros(size(vel_zshift,1),abs(xshift),size(vel_zshift,3));
                vel_xshift = cat(2,vel_zshift(:,abs(xshift-1):end,:),xadd);
            else
                vel_xshift = vel_zshift;
            end
            % Complete y-shift
            % Adjust rows based on y shift value
            if yshift > 0
                yadd = zeros(yshift,size(vel_xshift,2),size(vel_xshift,3));
                vel_shifted = cat(1,yadd,vel_xshift(1:end-yshift,:,:));
            elseif yshift < 0
                yadd = zeros(abs(yshift),size(vel_xshift,2),size(vel_xshift,3));
                vel_shifted = cat(1,vel_xshift(abs(yshift-1):end,:,:),yadd);
            else
                vel_shifted = vel_xshift;
            end
            cadat.DATA.Ur = vel_shifted.*cadat.STL.maskf;
            cadat.DATA.Vr = vel_shifted.*cadat.STL.maskf;
            cadat.DATA.Wr = vel_shifted.*cadat.STL.maskf;
            [cadat.DATA.Xr,cadat.DATA.Yr,cadat.DATA.Zr] = meshgrid(stl_x,stl_y,stl_z);
        else
            % Each point needs to be determined if it is in or out of the
            % mask.
            % Convert the x, y, z coordinates based on the shift,
            % centering, and conversion to stl coordinates
            % Center coordinates
            xc = cadat.DATA.x - cadat.DATA.centvalX;
            yc = cadat.DATA.y - cadat.DATA.centvalY;
            zc = cadat.DATA.z - cadat.DATA.centvalZ;
            
            % Flip coordinates if the data was flipped> No! Should be
            % flipped already
            
            
            % Adjust based on shift
            xf = xc + xshift*cadat.STL.dx;
            yf = yc + yshift*cadat.STL.dy;
            zf = zc + zshift*cadat.STL.dz;
            
            % Adjust based on rotation
            xf = xf - cadat.DATA.medX;
            yf = yf - cadat.DATA.medY;
            zf = zf - cadat.DATA.medZ;
            cXYZ = [xf';yf';zf'];
            % X rotation
            % Get the X rotation matrix
            rotmat = rotx(cadat.DATA.rot_optX);
            % Administer the X rotation
            cXYZ_rotX = rotmat*cXYZ;
            % Get the X rotation matrix
            rotmat = roty(cadat.DATA.rot_optY);
            % Administer the Y rotation
            cXYZ_rotXY = rotmat*cXYZ_rotX;
            % Get the Z rotation matrix
            rotmat = rotz(cadat.DATA.rot_optZ);
            % Administer the Z rotation
            cXYZ_rotXYZ = rotmat*cXYZ_rotXY;
            
            % Add back the median value
            xf = transpose(cXYZ_rotXYZ(1,:)) + cadat.DATA.medX;
            yf = transpose(cXYZ_rotXYZ(2,:)) + cadat.DATA.medY;
            zf = transpose(cXYZ_rotXYZ(3,:)) + cadat.DATA.medZ;
            
            % For each point, determine if it is in the mask or out
            add_ct = 1;
            for zz = 1:1:length(cadat.DATA.x)
                % Get current x, y, z coordinates
                curr_x = xf(zz);
                curr_y = yf(zz);
                curr_z = zf(zz);
                
                % Determine if the current point is in the mask or not by
                % finding the nearest mask point. If that mask point is 1,
                % then the point stays, else it is removed.
                [~,r_ind] = min(abs(stl_y - curr_y));
                [~,c_ind] = min(abs(stl_x - curr_x));
                [~,z_ind] = min(abs(stl_z - curr_z));
                
                in_mask = cadat.STL.maskf(r_ind,c_ind,z_ind);
                if in_mask == 1
                    % Save the point for the output
                    cadat.DATA.Ur(add_ct) = cadat.DATA.u(zz);
                    cadat.DATA.Vr(add_ct) = cadat.DATA.v(zz);
                    cadat.DATA.Wr(add_ct) = cadat.DATA.w(zz);
                    cadat.DATA.Xr(add_ct) = curr_x;
                    cadat.DATA.Yr(add_ct) = curr_y;
                    cadat.DATA.Zr(add_ct) = curr_z;
                    add_ct = add_ct + 1;
                end
            end
        end
        fprintf('completed successfully') 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SAVE MASKED AND REGISTERED VELOCITY FIELD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Only form velocity lines if saving is to be done
    if (cadat.OPTIONS.SAVE_RM_DATS || cadat.OPTIONS.SAVE_MAT) && strcmp(cadat.FILE.filetype,'TOMO')
        % Reshape grid into lines
        cadat.DATA.Xr = reshape(cadat.DATA.Xr,[size(cadat.DATA.Xr,1)*size(cadat.DATA.Xr,2)*size(cadat.DATA.Xr,3),1,1]);
        cadat.DATA.Yr = reshape(cadat.DATA.Yr,[size(cadat.DATA.Yr,1)*size(cadat.DATA.Yr,2)*size(cadat.DATA.Yr,3),1,1]);
        cadat.DATA.Zr = reshape(cadat.DATA.Zr,[size(cadat.DATA.Zr,1)*size(cadat.DATA.Zr,2)*size(cadat.DATA.Zr,3),1,1]);

        % Reshape masked velocities into lines
        cadat.DATA.Ur = reshape(cadat.DATA.Ur,[size(cadat.DATA.Ur,1)*size(cadat.DATA.Ur,2)*size(cadat.DATA.Ur,3),1,1]);
        cadat.DATA.Vr = reshape(cadat.DATA.Vr,[size(cadat.DATA.Vr,1)*size(cadat.DATA.Vr,2)*size(cadat.DATA.Vr,3),1,1]);
        cadat.DATA.Wr = reshape(cadat.DATA.Wr,[size(cadat.DATA.Wr,1)*size(cadat.DATA.Wr,2)*size(cadat.DATA.Wr,3),1,1]);
    end
    %%% Save the .mat files
    if cadat.OPTIONS.SAVE_MAT
        fprintf('\nSaving MAT file for file %i...',vnum)
        % Get file name
        fname = [cadat.DIR.velout,cadat.FILE.outputfilebase,num2str(vnum,'%05i'),'.mat'];
        % Get output
        x = cadat.DATA.Xr;
        y = cadat.DATA.Yr;
        z = cadat.DATA.Zr;
        u = cadat.DATA.Ur;
        v = cadat.DATA.Vr;
        w = cadat.DATA.Wr;
        % Save to mat file
        save(fname,'x','y','z','u','v','w')
        fprintf('completed successfully') 
    end
    
    %%% Save the .dat files
    if cadat.OPTIONS.SAVE_RM_DATS
        fprintf('\nSaving DAT file for file %i...',vnum)
        % Initialize output file
        fid = fopen([cadat.DIR.velout,cadat.FILE.outputfilebase,num2str(vnum,'%05i'),'.dat'],'w');
        fprintf(fid,['TITLE = "',cadat.FILE.outputfilebase,num2str(vnum,'%05i'),'"']);
        fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "V-X", "V-Y", "V-Z"\n');
        if cadat.OPTIONS.SAVE_ZEROS_DAT % Print the grid size if the zeros are being saved
            fprintf(fid,['ZONE T="Frame 0", I=',num2str(size(cadat.DATA.Ur,1)),', J=',num2str(size(cadat.DATA.Ur,2)),', K=',num2str(size(cadat.DATA.Ur,3))'']);
        end

        % Iterate through all points and write only those that are non-zero
        for zz = 1:1:length(cadat.DATA.Ur)
            if (abs(cadat.DATA.Ur(zz)) > 0) || cadat.OPTIONS.SAVE_ZEROS_DAT
                cx = cadat.DATA.Xr(zz);
                cy = cadat.DATA.Yr(zz);
                cz = cadat.DATA.Zr(zz);
                cu = cadat.DATA.Ur(zz);
                cv = cadat.DATA.Vr(zz);
                cw = cadat.DATA.Wr(zz);
                fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
            end
        end
        % Close dat file
        fclose(fid);
        fprintf('completed successfully') 
    end
    % For the first and second time steps, plot the final results to the 
    % user to ensure everything is running correctly to prevent the whole
    % code form running incorrectly.
    if (vnum == cadat.TIME.ts) || (vnum == cadat.TIME.ts+1)
        stl_xplot = stl_x;
        stl_yplot = stl_y;
        stl_zplot = stl_z;
        
        % Plot the data mask with the STL mask
        % Plot the STL and each velocity field
        [Xs,Ys,Zs] = meshgrid(stl_x,stl_y,stl_z);
        stl_mask = cadat.STL.maskf;
        figure(31); close(gcf);
        figure(31);
        vv = double(stl_mask > 0);
        p2 = patch(isosurface(Xs,Ys,Zs,vv,0));
        p2.FaceColor = [0.7 0.7 0.7];
        p2.EdgeColor = 'none';
        p2.FaceAlpha = 0.2;
        view(-90,0);
        
        hold on; scatter3(cadat.DATA.Xr,cadat.DATA.Yr,cadat.DATA.Zr,4,'filled');
        view(-90,0); 
        set(gca,'FontSize',14); title('Registered Data');
        prompt = '\nRegistration is acceptable? [1-Yes, 0-No]: ';
        accept_registration = input(prompt);
        if accept_registration ~= 1
            fprintf('ERROR - REGISTRATION REJECTED!\n PLEASE CHECK AND ADJUST CODE\n')
            fprintf('PROGRAM BEING ABORTED!\n')
            return
        end
    end
    
    fprintf('\n---------- File %i completed successfuly ----------\n',vnum)
end

end % End of function

