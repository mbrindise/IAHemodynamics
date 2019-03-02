
%% DIRECTORY AND FILE NAME INFORMATION
cadat.FILE.anynum = '007'; % Aneurysm identifiable number
cadat.FILE.filetype = 'CFD'; % Should be STB (Shake the box), TOMO, MRI, or CFD
cadat.SYSTEM.cpu = 'pc'; % Can be mac, pc, or server - changes the directories appropriately
    % STB loads mat file, TOMO loads dat file
    % MRI loads dat file

%% SET OPTIONS AND CONSTANTS
% Experimental constants
cadat.DATA.rho = 1103; % In kg/m^3
cadat.DATA.mu = 3.5e-3; % In pa-s
cadat.DATA.nu = cadat.DATA.mu/cadat.DATA.rho; % In m^2/s

% Time step options for the current file type
cadat.options.TS = 1; % Start time
cadat.options.Tskip = 1; % Time skip (1 = every one, 2 = every other, etc.) (use 2 for ANY111 and 3 for ANY007)
cadat.options.TE = 387;%1321; % End time (ANY111 = 1112, ANY007 = 1321)
% Set the time step information for the MRI file type
cadat.options.MRI_TS = 1; % Start time for MRI
cadat.options.MRI_Tskip = 1; % Time skip for MRI
cadat.options.MRI_TE = 13; % End time for MRI (ANY111 = 20, ANY007 = 13)
% Set t_skip for Shake the Box processing
cadat.options.STB_Tskip = 3;
% Set dt for each time step
cadat.options.dt.mri = 0.0448; % For ANY111=.0405, for ANY007=0.0448
cadat.options.dt.stb = 1/2000 * cadat.options.STB_Tskip; % For ANY111=2*1/400, for ANY007=3*1/400, 
cadat.options.dt.cfd = 0.0015; % For ANY111=0.0015, for ANY007=0.0015

% Results and computation options
% Velocity smoothing and data set options
cadat.options.RUNPOD = 0; % Run POD on velocity fields (should use gridded, masked, and registered fields as input)
    cadat.options.USEPOD = 1; % Use the POD filtered velocity fields for all subsequent calculations
cadat.options.RUNUOD = 0; % Run UOD on velocity fields (should use POD filtered fields as input)
    cadat.options.USEUOD = 1; % Use the UOD filtered velocity fields for all subsequent calculations
cadat.options.RUNPHASEAVGSTB = 0; % Run phase averaging on the velocity fields (only needed for the PIV data)
    cadat.options.USEPHASEAVG = 1; % Use phase averaging for all subsequent calculations
cadat.options.RUNVOXAVG = 0; % Run "voxel averaging" for STB/CFD data for velocity comparisons to MRI
    cadat.options.USEVOXAVG = 0; % Use the voxel averaged velocity fields for all subsequent calculation
                                 % Voxel averaging should *ONLY* be used for velocity comparsion,
                                 % NOT wss, osi rrt, pressure, coherent structures, etc.
cadat.options.CONVERT3DGRID = 0; % Converts the grid to 3D for paraview plotting

% Velocity and flow rate comparisons and calculations
cadat.options.SAVE_VEL_STL_PLOT = 0; % Save the STL plot needed for plotting the WSS
cadat.options.COMPFLOWRATE = 0; % Compute the flow rate for the current modality
    cadat.options.SETXSECS = 1; % The cross sections can be loaded or computed
% Secondary flow variable (i.e. wss, pressure, etc.) calculations
cadat.options.COMPNORMS = 0; % Compute the wall normal vectors for WSS
    cadat.options.SMOOTHNORMS = 1; % Smooth the normals previously computed
    cadat.options.LOADNORMS = 1; % Load in pre-computed velocity normals
cadat.options.COMPVELGRAD = 0; % Compute velocity gradients
cadat.options.COMPDIVERGENCE = 0; % Compute divergence (requires gradients to have been computed already)
cadat.options.COMPWSS = 0; % Compute wall shear stress (requires gradients to have been computed already)
cadat.options.COMPOSI = 1; % Compute oscillatory shear index (requires WSS to have already been computed)
cadat.options.COMPRRT = 1; % Compute relative residence times (requires WSS,OSI to have already been computed)
cadat.options.COMPVORT = 0; % Compute vorticity (requires gradients to have been computed already)
cadat.options.IDCOHERSTRUC = 0; % Identify coherent structures
cadat.options.COMPPRESSURE = 0; % Compute pressure
cadat.options.COMPPWV = 0; % Compute the pressure wave velocity
% Save output options
cadat.options.SAVEOUTPUTS = 1; % Choose whether to save the outputs or only compute them

% Plotting options
cadat.options.PLOTSTREAMTRACES = 0;
cadat.options.PLOTVELCOMP = 0;
cadat.options.PLOTVELPROFILES = 0;
cadat.options.PLOTDIVERGENCE = 0;
cadat.options.PLOTWSS = 0;

% Debug plotting option
cadat.options.DEBUGPLOTS = 1;

%% Set the directory and file names based on user inputs
% Set the base directory name and slash based on the system being used
if strcmp(cadat.SYSTEM.cpu,'pc')
    cadat.DIR.dirbase = 'Z:\Projects\Cerebral_Aneurysm\ANY\';
    fslash = '\';
elseif strcmp(cadat.SYSTEM.cpu,'server')
    cadat.DIR.dirbase = '~/Projects/Cerebral_Aneurysm/ANY/';
    fslash = '/';
else
    cadat.DIR.dirbase = '~/ANY/';
    fslash = '/';
    addpath ~/Documents/Code/general_codes/multipurpose
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the use/run controls based on the chosen modality
if strcmp(cadat.FILE.filetype,'MRI')
    cadat.options.USEUOD = 0;
    cadat.options.USEPOD = 0;
    cadat.options.USEVOXAVG = 0;
    cadat.options.USEPHASEAVG = 0;
elseif strcmp(cadat.FILE.filetype,'CFD')
    cadat.options.USEUOD = 0;
    cadat.options.USEPOD = 0;
    cadat.options.USEPHASEAVG = 0;
end

if strcmp(cadat.FILE.filetype,'STB')
    cadat.DIR.velfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_gridded',fslash];
    cadat.FILE.basename = ['any',cadat.FILE.anynum,'_stb_vi_regmask_gridded_'];
    cadat.FILE.podname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_v1_p1_pod_'];
    cadat.FILE.uodname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_v1_p2_uod_'];
    cadat.FILE.paname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_v1_p2s1_phaseavg_'];
    cadat.FILE.voxavgname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_v1_p3_voxavg_'];
elseif strcmp(cadat.FILE.filetype,'MRI')
    %cadat.DIR.velfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_registered_masked',fslash];
    %cadat.FILE.basename = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_regmask_'];
    cadat.DIR.velfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_gridded',fslash];
    cadat.FILE.basename = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_regmask_gridded_'];
else
    cadat.DIR.velfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_gridded',fslash];
    cadat.FILE.basename = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_regmask_gridded_'];
    cadat.FILE.podname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_p1_pod_'];
    cadat.FILE.uodname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_p2_uod_'];
    cadat.FILE.voxavgname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_p3_voxavg_'];
end
cadat.FILE.grid3Dname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_vA_grid3D_'];
cadat.DIR.podfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_p1_pod',fslash];
cadat.DIR.uodfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_p2_uod',fslash];
cadat.DIR.pafiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_p2s1_phaseavg',fslash];
cadat.DIR.voxavgfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_p3_voxavg',fslash];
cadat.DIR.grid3Dfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_pA_grid3D',fslash];
cadat.DIR.stlfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,'geometry',fslash];
cadat.DIR.savefiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,'registration_files',fslash];
cadat.DIR.saveppfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,'postprocessing_files',fslash];
cadat.DIR.saveplots = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,'plots',fslash];
cadat.FILE.anyname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_'];
cadat.FILE.maskgrid = ['ANY-',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_masked_velocity_grid'];
cadat.FILE.stlmaskname = ['STL-ANY-',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_STLMask_filt'];


%% INITIAL CALCULATION
% Load the grid mask for the velocity field
% This step is needed for several subsequent processes
fp = load([cadat.DIR.savefiles,cadat.FILE.maskgrid,'.mat']);
velmask = fp.velmask; % Velocity grid mask
maskX = fp.xmask;
maskY = fp.ymask;
maskZ = fp.zmask;
[Vx,Vy,Vz] = meshgrid(maskX,maskY,maskZ);
dx = maskX(2) - maskX(1); % In mm
dy = maskY(2) - maskY(1); % In mm
dz = maskZ(2) - maskZ(1); % In mm
cadat.DATA.X = Vx;
cadat.DATA.Y = Vy;
cadat.DATA.Z = Vz;

% Load the STL mask  - needed for computing flow rate
fp = load([cadat.DIR.savefiles,cadat.FILE.stlmaskname,'.mat']);
cadat.STL.maskf = fp.stlmaskf;
cadat.STL.maskf(cadat.STL.maskf == 0) = NaN;
cadat.STL.x = fp.stl_x;
cadat.STL.y = fp.stl_y;
cadat.STL.z = fp.stl_z;
unqX = unique(cadat.STL.x);
cadat.STL.dx = unqX(2)-unqX(1);
unqY = unique(cadat.STL.y);
cadat.STL.dy = unqY(2)-unqY(1);
unqZ = unique(cadat.STL.z);
cadat.STL.dz = unqZ(2)-unqZ(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VELOCITY FILTERING OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PROPER ORTHOGONAL DECOMPOSITION SMOOTHING
% This step smooths the data using POD in 3D. If POD is to be used for the
% remainder of the code (based on USEPOD option), then it will replace the
% directory and file names of the velocity fields
if cadat.options.RUNPOD
    fprintf('\nRunning POD on data...\n')
    % Load in all velocity files, input needs to be a 4D matrix (x,y,z,t).
    % Note** This could mean a VERY BIG matrix that requires A LOT of RAM!
    % Be careful to ensure the processing computer has enough memory!
    %%% Load velocity file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Loading data...')
    iter = 1;
    prev_mark = 0;
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
        % Get current velocities
        x = fp.x;
        y = fp.y;
        z = fp.z;
        u = fp.u;
        v = fp.v;
        w = fp.w;
        t = fp.t;
        
        % Place the velocities on a grid
        cu = zeros(size(velmask));
        cv = zeros(size(velmask));
        cw = zeros(size(velmask));
        for ptn = 1:1:length(x)
            [~,cind] = min(abs(maskX - x(ptn)));
            [~,rind] = min(abs(maskY - y(ptn)));
            [~,zind] = min(abs(maskZ - z(ptn)));
            
            % Put velocity in corresponding spatial location
            cu(rind,cind,zind) = u(ptn);
            cv(rind,cind,zind) = v(ptn);
            cw(rind,cind,zind) = w(ptn);
        end
        % Save all velocities into single array
        ua(:,:,:,iter) = cu;
        va(:,:,:,iter) = cv;
        wa(:,:,:,iter) = cw;
        ta(iter,1) = t;
        % Increment count
        iter = iter + 1;
        
        % Print the percent completed for steps of 10
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
    end
    fprintf(' completed successfully')
    %%% Run POD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nPerforming POD...')
    [Us,Vs,Ws,mkp,D] = POD3(ua,va,wa,2);
    fprintf('\nPOD completed successfully')
    
    keyboard
    
    %%% Save POD files for each time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nSaving POD velocity files...')
    % Convert x, y, z grids to lines
    Vxl = reshape(Vx,[size(Vx,1)*size(Vx,2)*size(Vx,3),1]);
    Vyl = reshape(Vy,[size(Vx,1)*size(Vx,2)*size(Vx,3),1]);
    Vzl = reshape(Vz,[size(Vx,1)*size(Vx,2)*size(Vx,3),1]);
    
    % Iterate through all time steps to save to mat and dat files
    iter = 1;
    prev_mark = 0;
    for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE %cadat.options.TS
        % Get current velocity fields
        uc = Us(:,:,:,iter);
        vc = Vs(:,:,:,iter);
        wc = Ws(:,:,:,iter);
        
        % Mask to only include the velocity mask points
        uc( isnan(velmask) ) = NaN;
        vc( isnan(velmask) ) = NaN;
        wc( isnan(velmask) ) = NaN;
        
        % Convert velocities to line
        uc = reshape(uc,[size(uc,1)*size(uc,2)*size(uc,3),1]);
        vc = reshape(vc,[size(uc,1)*size(uc,2)*size(uc,3),1]);
        wc = reshape(wc,[size(uc,1)*size(uc,2)*size(uc,3),1]);
        
        % Keep only the real values (get rid of the NaNs)
        x = Vxl(~isnan(uc));
        y = Vyl(~isnan(uc));
        z = Vzl(~isnan(uc));
        u = uc(~isnan(uc));
        v = vc(~isnan(uc));
        w = wc(~isnan(uc));
        t = ta(iter,1);
        
        % Save mat file
        save([cadat.DIR.podfiles,cadat.FILE.podname,num2str(q,'%05i'),'.mat'],'x','y','z','u','v','w','t');
        
        % Save dat file
        % Initialize output file
        fid = fopen([cadat.DIR.podfiles,cadat.FILE.podname,num2str(q,'%05i'),'.dat'],'w');
        fprintf(fid,['TITLE = "',cadat.FILE.podname,num2str(q,'%05i'),'"']);
        fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "U", "V", "W"\n');

        % Iterate through all points and write only those that are non-zero
        for zz = 1:1:length(x)
            cx = x(zz);
            cy = y(zz);
            cz = z(zz);
            cu = u(zz);
            cv = v(zz);
            cw = w(zz);
            fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
        end
        % Close dat file
        fclose(fid);
        
        % Increment count
        iter = iter + 1;
        
        % Print the percent completed for steps of 10
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
        
    end
    fprintf('\n...completed successfully\n')
end

if cadat.options.USEPOD
    % Set directory and base name of velocity files to the pod ones
    cadat.DIR.velfiles = cadat.DIR.podfiles;
    cadat.FILE.basename = cadat.FILE.podname;
end

%% UNIVERSAL OUTLIER DETECTION
% This step runs a 3D UOD on each velocity field in time. If UOD is to be
% used for the remainder of the code (based on USEUOD option), then it will
% replace the directory and file names of the velocity field
if cadat.options.RUNUOD
    fprintf('\nRunning UOD on data...')
    uodmask = ~isnan(velmask);
    % Convert x, y, z grids to lines
    Vxl = reshape(Vx,[size(Vx,1)*size(Vx,2)*size(Vx,3),1]);
    Vyl = reshape(Vy,[size(Vx,1)*size(Vx,2)*size(Vx,3),1]);
    Vzl = reshape(Vz,[size(Vx,1)*size(Vx,2)*size(Vx,3),1]);
    
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    iter = 1;
    prev_mark = 0;
    format long
    for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        %%% Load velocity file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
        % Get current velocities
        x = fp.x;
        y = fp.y;
        z = fp.z;
        u = fp.u;
        v = fp.v;
        w = fp.w;
        t = fp.t;
        % Place the velocities on a grid
        cu = zeros(size(velmask));
        cv = zeros(size(velmask));
        cw = zeros(size(velmask));
        for ptn = 1:1:length(x)
            [~,cind] = min(abs(maskX - x(ptn)));
            [~,rind] = min(abs(maskY - y(ptn)));
            [~,zind] = min(abs(maskZ - z(ptn)));
            
            % Put velocity in corresponding spatial location
            cu(rind,cind,zind) = u(ptn);
            cv(rind,cind,zind) = v(ptn);
            cw(rind,cind,zind) = w(ptn);
        end
        
        %%% Run UOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [uf,vf,wf] = UOD3(cu,cv,cw,uodmask,5,2.5,50,dx);
        
        %%% Save UOD mat and dat files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get current velocity fields
        uc = uf;%(:,:,:,iter);
        vc = vf;%(:,:,:,iter);
        wc = wf;%(:,:,:,iter);

        % Mask to only include the velocity mask points
        uc( isnan(velmask) ) = NaN;
        vc( isnan(velmask) ) = NaN;
        wc( isnan(velmask) ) = NaN;

        % Convert velocities to line
        uc = reshape(uc,[size(uc,1)*size(uc,2)*size(uc,3),1]);
        vc = reshape(vc,[size(uc,1)*size(uc,2)*size(uc,3),1]);
        wc = reshape(wc,[size(uc,1)*size(uc,2)*size(uc,3),1]);

        % Keep only the real values (get rid of the NaNs)
        x = Vxl(~isnan(uc));
        y = Vyl(~isnan(uc));
        z = Vzl(~isnan(uc));
        u = uc(~isnan(uc));
        v = vc(~isnan(uc));
        w = wc(~isnan(uc));

        % Save mat file
        save([cadat.DIR.uodfiles,cadat.FILE.uodname,num2str(q,'%05i'),'.mat'],'x','y','z','u','v','w','t');

        % Save dat file
        % Initialize output file
        fid = fopen([cadat.DIR.uodfiles,cadat.FILE.uodname,num2str(q,'%05i'),'.dat'],'w');
        fprintf(fid,['TITLE = "',cadat.FILE.uodname,num2str(q,'%05i'),'"']);
        fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "U", "V", "W"\n');

        % Iterate through all points and write only those that are non-zero
        for zz = 1:1:length(x)
            cx = x(zz);
            cy = y(zz);
            cz = z(zz);
            cu = u(zz);
            cv = v(zz);
            cw = w(zz);
            fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
        end
        % Close dat file
        fclose(fid);

        % Increment count
        iter = iter + 1;

        % Print the percent completed for steps of 5
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,5) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
    end
    fprintf(' completed successfully')
    format short
end

if cadat.options.USEUOD
    % Set directory and base name of velocity files to the pod ones
    cadat.DIR.velfiles = cadat.DIR.uodfiles;
    cadat.FILE.basename = cadat.FILE.uodname;
end


%% PHASE AVERAGE PIV DATA
if cadat.options.RUNPHASEAVGSTB
    fprintf('\nRunning phase averaging on data...\n')
    %%% Load velocity file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q = cadat.options.TS;
    fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
    
    %%% Have the user select a point to use for the phase averaging %%%%%%%
    % This velocity at this point is extracted
    figure(51); hold off; scatter3(fp.x,fp.y,fp.z); view(2); set(gcf,'Color',[1 1 1]);
    set(gca,'FontSize',14); title('Scatter of grid points'); xlabel('X'); ylabel('Y');
    % Get the X, Y, Z position where the point velocity should be
    prompt = 'Enter [X,Y,Z] for point velocity extraction: ';
    pt_xsec = input(prompt);
    % Save the points for each cross section
    cx = pt_xsec(1);
    cy = pt_xsec(2);
    cz = pt_xsec(3);

    % Find the x and y indices of the current point
    [~,cxi] = min(abs(maskX - cx));
    [~,cyi] = min(abs(maskY - cy));
    [~,czi] = min(abs(maskZ - cz));
    % Find the STL x and y indices of the current point
    [~,cxSTLi] = min(abs(cadat.STL.x - cx));
    [~,cySTLi] = min(abs(cadat.STL.y - cy));
    [~,czSTLi] = min(abs(cadat.STL.z - cz));

    % Get the z line corresponding to the x and y points and find
    % where the mask is defined
    z_line = reshape(cadat.STL.maskf(cySTLi,cxSTLi,:),[size(cadat.STL.maskf,3),1]);
    zs = find(isnan(z_line(1:czSTLi)),1,'last') + 1;
    ze = find(isnan(z_line(czSTLi:end)),1,'first') + czSTLi - 2;
    % Ensure zs and ze are defined
    if isempty(zs), zs = 1; end
    if isempty(ze), ze = length(x_line); end
    % Get the x line corresponding to the z and y points and find
    % where the mask is defined
    x_line = reshape(cadat.STL.maskf(cySTLi,:,czSTLi),[size(cadat.STL.maskf,2),1]);
    xs = find(isnan(x_line(1:cxSTLi)),1,'last') + 1;
    xe = find(isnan(x_line(cxSTLi:end)),1,'first') + cxSTLi - 2;
    % Ensure xs and xe are defined
    if isempty(xs), xs = 1; end
    if isempty(xe), xe = length(x_line); end
    % Get the y line corresponding to the x and z points and find
    % where the mask is defined
    y_line = reshape(cadat.STL.maskf(:,cxSTLi,czSTLi),[size(cadat.STL.maskf,1),1]);
    ys = find(isnan(y_line(1:cySTLi)),1,'last') + 1;
    ye = find(isnan(y_line(cySTLi:end)),1,'first') + cySTLi - 2;
    % Ensure ys and ye are defined
    if isempty(ys), ys = 1; end
    if isempty(ye), ye = length(y_line); end
    % Get x and y and z distances
    x_dist = xe - xs;
    y_dist = ye - ys;
    z_dist = ze - zs;
    
    % Find the centroid of the x, y, z point based on the distances
    % computed
    if (x_dist > y_dist) && (x_dist > z_dist)
        xcent = cxSTLi;
        ycent = round(ys + (ye-ys)/2);
        zcent = round(zs + (ze-zs)/2);
    elseif (y_dist > x_dist) && (y_dist > z_dist)
        xcent = round(xs + (xe-xs)/2);
        ycent = cySTLi;
        zcent = round(zs + (ze-zs)/2);
    else
        xcent = round(xs + (xe-xs)/2);
        ycent = round(ys + (ye-ys)/2);
        zcent = czSTLi;
    end
    
    % Get the physical coordinates of the chosen point
    x_vpt = cadat.STL.x(xcent);
    y_vpt = cadat.STL.y(ycent);
    z_vpt = cadat.STL.z(zcent);
    
    % Iterate through all time and collect the median velocity of the 5
    % points closest to the center point of the cross section identified by
    % the user
    % Initialize velocity array
%     u_vpt = zeros(length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE),1);
%     v_vpt = zeros(length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE),1);
%     w_vpt = zeros(length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE),1);
%     % Find the 5 points closest to the point with physical coordinates
%     ct = 1;
%     for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
%         %%% Load velocity file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
%         % Get current velocities
%         x = fp.x; y = fp.y; z = fp.z;
%         u = fp.u; v = fp.v; w = fp.w;
%         t = fp.t;
%         % Compute the distance of each point from the selected point
%         pt_dists = sqrt((x-x_vpt).^2 + (y-y_vpt).^2 + (z-z_vpt).^2);
%         % Sort the points based on the distances
%         [srt_dist,srt_pts] = sort(pt_dists,'ascend');
%         % Add the point to the array
%         u_vpt(ct) = median(u(srt_pts(1:5)));
%         v_vpt(ct) = median(v(srt_pts(1:5)));
%         w_vpt(ct) = median(w(srt_pts(1:5)));
%         % Increment the counter
%         ct = ct + 1;
%     end
%     % Compute the velocity magnitude of the point
%     vel_vpt_raw = sqrt(u_vpt.^2 + v_vpt.^2 + w_vpt.^2);
%     T = wpdec(vel_vpt_raw,5,'sym4');
%     vel_vpt = wprcoef(T,[5 0]);
    
    %%% Use the velocity of the identified point to find the best overlap
    %%% of the cycles
    % Have the user input the number of complete cycles
    figure(52); plot(vel_vpt,'LineWidth',1.5);
    prompt = 'Enter the number of complete cycles: ';
    num_cyc = input(prompt);
    prompt = 'Enter the start index of the first complete cycle: ';
    cyc1S = input(prompt);
    prompt = 'Enter the end index of the first complete cycle: ';
    cyc1E = input(prompt);
    % Obtain the velocity of the first cycle
    cyc1_vel = vel_vpt(cyc1S:cyc1E);
    mid_cyc1 = round((cyc1S+cyc1E)/2);
    rm_vel = vel_vpt(mid_cyc1:end);
    [r,lags] = xcorr(vel_vpt(mid_cyc1:end),cyc1_vel);
    
    % Find the correlation peaks
    [pks,locs,wds,proms] = findpeaks(r); % Compute peaks of the correlation vector
    [pks,locs,wds,proms] = findpeaks(r,'MinPeakDistance',(cyc1E-cyc1S)*0.75); % Compute peaks of the correlation vector
    num_pks = num_cyc-1; % Number of peaks to correlate, based on the number of cycles
    % Sort the peaks based on the peak prominence
    [srt_prom,srt_pts] = sort(proms,'descend');
    cyc_locs = locs(srt_pts(1:num_pks));
    % From the locations, the starting points of the remaining cycles can
    % be computed from the cycle locations
    cyc_lags = sort(lags(cyc_locs),'ascend'); % Compute the lags for the remaining cycle starts
    %cyc_lags = [175,526,912];
    cyc_spts = cyc_lags + mid_cyc1; % Compute the cycle start points for the remaining cycles
    cycS = [cyc1S,cyc_spts]; % Obtain the cycle start points for each cycle
    % Compute the end of the last cycle
    cyc_len = cycS(2:end)-cycS(1:end-1);
    cyc_med_len = round(median(cyc_len));
    cycE = [cycS(2:end)-1,cycS(end)+cyc_med_len]; % Array containing the end of the cycles
    
    % Use the cycle length as the median of the cycles
    cyc_len = round(median(cycE - cycS));
    
    % Obtain the time points of each cycle
    cyc_time_pts = zeros(num_cyc,cyc_len);
    for ncyc = 1:1:num_cyc
        cyc_pts = round(linspace(cycS(ncyc),cycE(ncyc),cyc_len));
        cyc_time_pts(ncyc,:) = cyc_pts;
    end
    
    % Convert the time points to the file numbers
    file_nums = cadat.options.TS:cadat.options.Tskip:cadat.options.TE;
    cyc_file_pts = file_nums(cyc_time_pts);
    
    % Plot the phase averaging example using the velocity
    figure(53); hold off; plot(vel_vpt(cycS(1):cycE(1)),'k','LineWidth',1.5); hold on;
    for nn = 2:num_cyc
        plot(vel_vpt(cycS(nn):cycE(nn)),'k','LineWidth',1.5);
    end
    vel_avg = vel_vpt(cyc_time_pts);
    vel_avg = mean(vel_avg,1);
    t_cvel = 0:cadat.options.dt.stb:cadat.options.dt.stb*(length(vel_avg)-1);
    plot(vel_avg,'r','LineWidth',1.5);
    title('Phase Averaged Velocity Point');
    fprintf('...phase averaged cycle files computed successfully')
    
    %%% Determine the re-arrangement of the data so the cycle matches the
    %%% MRI cycle (CFD will match the MRI already) %%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize velocity array
    u_vpt = zeros(length(cadat.options.MRI_TS:cadat.options.MRI_Tskip:cadat.options.MRI_TE),1);
    v_vpt = zeros(length(cadat.options.MRI_TS:cadat.options.MRI_Tskip:cadat.options.MRI_TE),1);
    w_vpt = zeros(length(cadat.options.MRI_TS:cadat.options.MRI_Tskip:cadat.options.MRI_TE),1);
    % Find the 5 points closest to the point with physical coordinates
    ct = 1;
    % Load the first MRI velocity field
    mri_velfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,'MRI',fslash,'vel_gridded',fslash];
    mri_basename = ['any',cadat.FILE.anynum,'_','MRI','_regmask_gridded_'];
    for q = cadat.options.MRI_TS:cadat.options.MRI_Tskip:cadat.options.MRI_TE
        fp = load([mri_velfiles,mri_basename,num2str(q,'%05i'),'.mat']);
        % Get current velocities
        x = fp.x; y = fp.y; z = fp.z;
        u = fp.u; v = fp.v; w = fp.w;
        % Extract the velocites at the same point selected for the PIV data
        % Compute the distance of each point from the selected point
        pt_dists = sqrt((x-x_vpt).^2 + (y-y_vpt).^2 + (z-z_vpt).^2);
        % Sort the points based on the distances
        [srt_dist,srt_pts] = sort(pt_dists,'ascend');
        % Add the point to the array
        u_vpt(ct) = median(u(srt_pts(1:3)));
        v_vpt(ct) = median(v(srt_pts(1:3)));
        w_vpt(ct) = median(w(srt_pts(1:3)));
        % Increment the counter
        ct = ct + 1;
    end
    % Compute the velocity magnitude of the point
    vel_vpt = sqrt(u_vpt.^2 + v_vpt.^2 + w_vpt.^2);
    t_mri = 0:cadat.options.dt.mri:cadat.options.dt.mri*(length(vel_vpt)-1);
    
    % Interpolate the MRI velocity point to the same time resolution as the
    % current data set
    vel_int = interp1(t_mri,vel_vpt,t_cvel,'spline');
    vel_mri = vel_int;%[vel_int,vel_int];
    vel_mri = vel_mri/max(vel_mri);
    % Deteremine the optimal shift of the phase averaged data
    % Add an exra half cycle to the end of the PIV data
    vel_avg_add = [vel_avg,vel_avg];
    vel_avg_add = vel_avg_add/max(vel_avg_add);
    vel_avg_add_sm = smooth(vel_avg_add,'moving',7);
    [r,lags] = xcorr(vel_avg_add_sm,vel_mri);
    % Take only the positive lag values
    rc = r(lags > 0);
    lagsc = lags(lags > 0);
    % Find the maximum peak of the correlation
    [pks,locs,wds,proms] = findpeaks(rc);
    [~,optpk] = max(proms);
    opt_lag = locs(optpk);
    
    % Take the correct portion of the cycle time points based on the
    % optimal lag
    cyc_file_pts_rep = repmat(cyc_file_pts,[1,3]);
    cyc_time_pts_rep = repmat(cyc_time_pts,[1,3]);
    cyc_file_pts = cyc_file_pts_rep(:,opt_lag:opt_lag+cyc_len-1);
    cyc_time_pts = cyc_time_pts_rep(:,opt_lag:opt_lag+cyc_len-1);
    
    % Plot the correct orientation of the phase for the user
    vel_avg_add = [vel_avg_add,vel_avg];
    figure(54); hold off; plot(vel_mri,'k','LineWidth',1.5);
    hold on; plot(vel_avg_add(opt_lag:opt_lag+cyc_len-1),'r','LineWidth',1.5);
    
    %%% Complete phase averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nPhase averaging and saving data...')
    % Iterate through all phase averaged fields and complete phase
    % averaging
    clear XYZ UVW Time
    iter = 1; num_iter = size(cyc_file_pts,2); prev_mark = 0;
    for q = 1:1:size(cyc_file_pts,2)
        % Get the current file number (want to keep the same time skip as
        % before)
        curr_fp = q;%file_nums(q);
        
        % Load the current velocity fields to phase average  
        format long
        for nn = 1:1:num_cyc
            % Current file number to phase average
            curr_file = cyc_file_pts(nn,q);
            fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(curr_file,'%05i'),'.mat']);
            % Get current point locations and velocities
            x = fp.x; y = fp.y; z = fp.z;
            u = fp.u; v = fp.v; w = fp.w;
            curr_t = fp.t;
            % Add the point to the array
            XYZ{nn,1} = [x,y,z];
            UVW{nn,1} = [u,v,w];
        end
        
        % Rearrange the velocity fields into X,Y,Z arrays
        X = []; Y = []; Z = []; U = []; V = []; W = [];
        for nn = 1:1:num_cyc
            X = [X,XYZ{nn,1}(:,1)];
            Y = [Y,XYZ{nn,1}(:,2)];
            Z = [Z,XYZ{nn,1}(:,3)];
            U = [U,UVW{nn,1}(:,1)];
            V = [V,UVW{nn,1}(:,2)];
            W = [W,UVW{nn,1}(:,3)];
        end
        % Compute the time for the current file
        curr_time = (iter-1)*cadat.options.dt.stb;
        
        % Phase average the velocity fields
        Xpa = mean(X,2);
        Ypa = mean(Y,2);
        Zpa = mean(Z,2);
        Upa = mean(U,2);
        Vpa = mean(V,2);
        Wpa = mean(W,2);
        % Ensure that the X,Y,Z points align
        x_err = X - repmat(Xpa,[1,num_cyc]); x_err = sum(abs(x_err(:)));
        y_err = Y - repmat(Ypa,[1,num_cyc]); y_err = sum(abs(y_err(:)));
        z_err = Z - repmat(Zpa,[1,num_cyc]); z_err = sum(abs(z_err(:)));
        if (x_err + y_err + z_err) > 1e-6
            keyboard
        end
        
        %%% Export the phase averaged velocity fields %%%%%%%%%%%%%%%%%%%%%
        x = Xpa; y = Ypa; z = Zpa;
        u = Upa; v = Vpa; w = Wpa;
        t = curr_time;
        
        format short
        % Save mat file
        save([cadat.DIR.pafiles,cadat.FILE.paname,num2str(curr_fp,'%05i'),'.mat'],'x','y','z','u','v','w','t');

        % Save dat file
        if curr_fp <= 5
            % Initialize output file
            fid = fopen([cadat.DIR.pafiles,cadat.FILE.paname,num2str(curr_fp,'%05i'),'.dat'],'w');
            fprintf(fid,['TITLE = "',cadat.FILE.uodname,num2str(curr_fp,'%05i'),'"']);
            fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "U", "V", "W"\n');

            % Iterate through all points and write only those that are non-zero
            for zz = 1:1:length(x)
                cx = x(zz);
                cy = y(zz);
                cz = z(zz);
                cu = u(zz);
                cv = v(zz);
                cw = w(zz);
                fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
            end
            % Close dat file
            fclose(fid);

            % Print the percent completed for steps of 10
            perc_complete = floor(100*iter/num_iter);
            if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                fprintf(' %i%%',perc_complete)
                prev_mark = perc_complete;
            end
        end
        % Increment the iteration number
        iter = iter + 1;
    end
    % Save the number of phase averaged file outputs
    phase_avg_files = file_nums(1:iter-1);
    save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'phase_avg_info','.mat'],'phase_avg_files','cyc_file_pts','cyc_time_pts','pt_xsec','num_cyc','cyc1S','cyc1E')
    fprintf(' completed successfully\n')
    format short
end
if cadat.options.USEPHASEAVG
    % Set directory and base name of velocity files to the pod ones
    cadat.DIR.velfiles = cadat.DIR.pafiles;
    cadat.FILE.basename = cadat.FILE.paname;
    % Adjust the time options for the current time step based on the phase
    % average files
    fp = load([cadat.DIR.saveppfiles,cadat.FILE.anyname,'phase_avg_info','.mat']);
    phase_avg_files = fp.phase_avg_files;
    cadat.options.TS = 1; %phase_avg_files(1);
    cadat.options.TE = length(phase_avg_files);%(end);
    cadat.options.Tskip = 1; %phase_avg_files(2)-phase_avg_files(1);
end


%% VOXEL AVERAGING
% This step runs a voxel averaging code to best compare the MRI and CFD/STB
% velocity fields and flow rates
if cadat.options.RUNVOXAVG
    fprintf('\nRunning voxel averaging on data...')
    % Load the first MRI velocity field
    mri_velfiles = [cadat.DIR.dirbase,cadat.FILE.anynum,fslash,'MRI',fslash,'vel_registered_masked',fslash];
    mri_basename = ['any',cadat.FILE.anynum,'_','MRI','_regmask_'];
    % Load the first MRI file
    fpm = load([mri_velfiles,mri_basename,num2str(1,'%05i'),'.mat']);
    mri_x = fpm.x;
    mri_y = fpm.y;
    mri_z = fpm.z;
    unq1 = unique(mri_x);
    dx = unq1(2)-unq1(1);
    unq1 = unique(mri_y);
    dy = unq1(2)-unq1(1);
    unq1 = unique(mri_z);
    dz = unq1(2)-unq1(1);
    
    % Compute the number of time steps to compute
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    iter = 1; prev_mark = 0;
    % Iterate through all time steps
    for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        % Load the current velocity file
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
        % Compute the x, y, z arrays for the current velocity file that
        % have the desired MRI dx, dy, dz
        x = fp.x; y = fp.y; z = fp.z;
        u = fp.u; v = fp.v; w = fp.w;
        t = fp.t;
        % Create the mesh grid for the voxel averaging - only do this for
        % the first time step
        if q == cadat.options.TS
            xL = min(x):dx:max(x);
            yL = min(y):dy:max(y);
            zL = min(z):dz:max(z);
            [Xmri,Ymri,Zmri] = meshgrid(xL,yL,zL);
            xm = reshape(Xmri,[size(Xmri,1)*size(Xmri,2)*size(Xmri,3),1]);
            ym = reshape(Ymri,[size(Xmri,1)*size(Xmri,2)*size(Xmri,3),1]);
            zm = reshape(Zmri,[size(Xmri,1)*size(Xmri,2)*size(Xmri,3),1]);
            velmask_va = nan(size(xm));
            % Determine which points are in the mask, keep only those points
            in_mask = zeros(size(xm));
            for p = 1:1:length(xm)
                % Obtain x, y, z of current grid points
                cx = xm(p);
                cy = ym(p);
                cz = zm(p);
                % Determine if the current point is in the velocity mask
                [~,c_ind] = min(abs(maskX-cx));
                [~,r_ind] = min(abs(maskY-cy));
                [~,z_ind] = min(abs(maskZ-cz));
                is_inmask = ~isnan(velmask(r_ind,c_ind,z_ind));
                in_mask(p) = is_inmask;
            end
            % Truncate to only the points in the mask
            xm_mask = xm(in_mask == 1);
            ym_mask = ym(in_mask == 1);
            zm_mask = zm(in_mask == 1);
            velmask_va(in_mask == 1) = 0;
            % Reshape the velocity mask
            velmask_va = reshape(velmask_va,[size(Xmri,1),size(Xmri,2),size(Xmri,3)]);
        end
        
        % Assign each STB/CFD point to the new grid point it is closest to
        pt_assign = cell(size(xm_mask));
        for p = 1:1:length(x)
            % Get the current STB/CFD point
            cx = x(p);
            cy = y(p);
            cz = z(p);
            % Compute its distance from each grid point in the mask
            pt_dist = sqrt((xm_mask - cx).^2 + (ym_mask - cy).^2 + (zm_mask - cz).^2);
            [~,min_dist] = min(pt_dist);
            % Obtain the cell of the closest index
            pt_cell = pt_assign{min_dist};
            % Append the current point to the cell
            pt_cella = [pt_cell;p];
            % Reset the cell to the appended matrix
            pt_assign{min_dist} = pt_cella;
        end
        
        % Iterate through each grid point and compute the velocities based
        % on a inverse radius weighted average
        % Initialize the output
        cu = zeros(size(xm_mask));
        cv = zeros(size(xm_mask));
        cw = zeros(size(xm_mask));
        for p = 1:1:length(xm_mask)
            % Get the current grid point
            cx = xm_mask(p);
            cy = ym_mask(p);
            cz = zm_mask(p);
            % Get the indices assigned to this grid point
            vel_inds = pt_assign{p};
            % Get the x,y,z,u,v,w coordinates assigned to this grid point
            in_rad_x = x(vel_inds);
            in_rad_y = y(vel_inds);
            in_rad_z = z(vel_inds);
            in_rad_u = u(vel_inds);
            in_rad_v = v(vel_inds);
            in_rad_w = w(vel_inds);
            % Compute radial distance for all points assigned to this grid
            % point
            in_rad_dists = sqrt((in_rad_x - cx).^2 + (in_rad_y - cy).^2 + (in_rad_z - cz).^2);

            % Compute weight based on square of distance
            dist_wgt = 1./(in_rad_dists.^2);
            dist_wgt = dist_wgt/sum(dist_wgt);

            % Save weighted velocity for point
            cu(p) = sum(in_rad_u.*dist_wgt);
            cv(p) = sum(in_rad_v.*dist_wgt);
            cw(p) = sum(in_rad_w.*dist_wgt);
        end
        % Get the output vectors
        x = xm_mask;
        y = ym_mask;
        z = zm_mask;
        u = cu;
        v = cv;
        w = cw;
        
        % Save mat file
        save([cadat.DIR.voxavgfiles,cadat.FILE.voxavgname,num2str(q,'%05i'),'.mat'],'x','y','z','u','v','w','t');

        % Save dat file
        % Initialize output file
        fid = fopen([cadat.DIR.voxavgfiles,cadat.FILE.voxavgname,num2str(q,'%05i'),'.dat'],'w');
        fprintf(fid,['TITLE = "',cadat.FILE.voxavgname,num2str(q,'%05i'),'"']);
        fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "U", "V", "W"\n');

        % Iterate through all points and write only those that are non-zero
        for zz = 1:1:length(x)
            cx = x(zz);
            cy = y(zz);
            cz = z(zz);
            cu = u(zz);
            cv = v(zz);
            cw = w(zz);
            fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
        end
        % Close dat file
        fclose(fid);

        % Increment count
        iter = iter + 1;

        % Print the percent completed for steps of 5
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,5) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
    end
    % Update the velocity mask and save the voxel averaged velocity mask
    velmask = velmask_va;
    xmask = xL; ymask = yL; zmask = zL;
    maskX = xL; maskY = yL; maskZ = zL;
    save([cadat.DIR.savefiles,cadat.FILE.maskgrid,'-voxel_averaged.mat'],'velmask','xmask','ymask','zmask')
    
    fprintf(' complete\n')
end
if cadat.options.USEVOXAVG
    % Set directory and base name of velocity files to the pod ones
    cadat.DIR.velfiles = cadat.DIR.voxavgfiles;
    cadat.FILE.basename = cadat.FILE.voxavgname;
    % Load the new velocity mask
    fp = load([cadat.DIR.savefiles,cadat.FILE.maskgrid,'-voxel_averaged.mat']);
    velmask = fp.velmask;
    maskX = fp.xmask;
    maskY = fp.ymask;
    maskZ = fp.zmask;
    [Vx,Vy,Vz] = meshgrid(maskX,maskY,maskZ);
    dx = maskX(2) - maskX(1); % In mm
    dy = maskY(2) - maskY(1); % In mm
    dz = maskZ(2) - maskZ(1); % In mm
end

%% CONVERT TO 3D GRID
% This is needed for ParaView where a complete 3D grid is needed
% Iterate through all time
if cadat.options.CONVERT3DGRID
    prev_mark = 0;
    iter = 1;
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    fprintf('\nSaving 3D grid data to dat files...')
    for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        %%% Load velocity file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
        % Form the grid
        xl = unique(fp.x);
        yl = unique(fp.y);
        zl = unique(fp.z);
        [Xg,Yg,Zg] = meshgrid(xl,yl,zl);
        I = size(Xg,1); J = size(Xg,2); K = size(Xg,3);
        Ug = zeros(size(Xg));
        Vg = zeros(size(Yg));
        Wg = zeros(size(Zg));
        % Place the velocities on the grid
        for zz = 1:1:length(fp.u)
            [~,c_ind] = min(abs(xl-fp.x(zz)));
            [~,r_ind] = min(abs(yl-fp.y(zz)));
            [~,z_ind] = min(abs(zl-fp.z(zz)));
            Ug(r_ind,c_ind,z_ind) = fp.u(zz);
            Vg(r_ind,c_ind,z_ind) = fp.v(zz);
            Wg(r_ind,c_ind,z_ind) = fp.w(zz);
        end
        % Reshape the arrays to 1D vectors
        x = reshape(Xg,[size(Xg,1)*size(Xg,2)*size(Xg,3),1]);
        y = reshape(Yg,[size(Xg,1)*size(Xg,2)*size(Xg,3),1]);
        z = reshape(Zg,[size(Xg,1)*size(Xg,2)*size(Xg,3),1]);
        u = reshape(Ug,[size(Xg,1)*size(Xg,2)*size(Xg,3),1]);
        v = reshape(Vg,[size(Xg,1)*size(Xg,2)*size(Xg,3),1]);
        w = reshape(Wg,[size(Xg,1)*size(Xg,2)*size(Xg,3),1]);
        % Save the 3D grid in a dat file
        % Save dat file
        % Initialize output file
        fid = fopen([cadat.DIR.grid3Dfiles,cadat.FILE.grid3Dname,num2str(q,'%05i'),'.dat'],'w');
        fprintf(fid,['TITLE = "',cadat.FILE.grid3Dname,num2str(q,'%05i'),'"']);
        fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "V-X", "V-Y", "V-Z"');
        fprintf(fid,['\nZONE T="FRAME 0", I=',num2str(I),', J=',num2str(J),', K=',num2str(K)]);
        % Iterate through all points and write only those that are non-zero
        for zz = 1:1:length(x)
            cx = x(zz);
            cy = y(zz);
            cz = z(zz);
            cu = u(zz);
            cv = v(zz);
            cw = w(zz);
            fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
        end
        % Close dat file
        fclose(fid);

        % Increment count
        iter = iter + 1;

        % Print the percent completed for steps of 5
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,5) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end   
    end
    fprintf(' complete\n')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VELOCITY COMPARISON OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE STREAMTRACES
% Computes streamtraces for plotting purposes for the current modality and
% velocity field selection
if cadat.options.SAVE_VEL_STL_PLOT
    save_vel_stl_plot(cadat);
end


% Change the name of the files based on full vs. voxel averaged resolution
if (strcmp(cadat.FILE.filetype,'CFD') || strcmp(cadat.FILE.filetype,'STB')) && cadat.options.USEVOXAVG == 0
    cadat.FILE.anyname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_full_resolution_'];
end

%% COMPUTE FLOW RATE
% This step computes the flow rate at specified cross sectional points. The
% cross sectional points can be loaded or computed as specified by the
% suboption.
% The cross sectional areas are selected using a frontal view 
if cadat.options.COMPFLOWRATE
    % Find the location of the flow rate estimates desired or import the
    % points if they have been previously defined
    if cadat.options.SETXSECS
        clear xsecX xsecY cross_sections 
        % Load in a single velocity field (the first one) and create a
        % single planed velocity mask
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(cadat.options.TS,'%05i'),'.mat']);
        % Plot data to user to receive input
        figure(51); hold off; scatter3(fp.x,fp.y,fp.z,4,'filled'); view(2); set(gcf,'Color',[1 1 1]);
        set(gca,'FontSize',14); title('Scatter of grid points'); xlabel('X'); ylabel('Y');
        % Get the number cross sections to be evaluated
        prompt = '\nHow many cross sections need to be defined?: ';
        num_xsecs = input(prompt);
        cadat.DATA.cross_section = cell(num_xsecs,1);
        % Iterate through all cross sections, for each have the user select
        % a point in the center of the vessel where the cross section is to
        % be computed
        for q = 1:1:num_xsecs
            % Get the X and Y positions where the cross section should be
            prompt = ['Enter [X,Y,Z] for cross section ',num2str(q),': '];
            pt_xsec = input(prompt);
            % Save the points for each cross section
            xsecX(q,1) = pt_xsec(1);
            xsecY(q,1) = pt_xsec(2);
            xsecZ(q,1) = pt_xsec(3);
        end
        
        % For each cross sectional point, need to determine the line that
        % goes perpendicularly across X,Y. Then find the angle
        % (z-direction)
        for q = 1:1:num_xsecs
            % Get x and y of current cross section center
            cx = xsecX(q);
            cy = xsecY(q);
            cz = xsecZ(q);
            
            % Find the x and y indices of the current point
            [~,cxi] = min(abs(maskX - cx));
            [~,cyi] = min(abs(maskY - cy));
            [~,czi] = min(abs(maskZ - cz));
            % Find the STL x and y indices of the current point
            [~,cxSTLi] = min(abs(cadat.STL.x - cx));
            [~,cySTLi] = min(abs(cadat.STL.y - cy));
            [~,czSTLi] = min(abs(cadat.STL.z - cz));
            
            % Get the z line corresponding to the x and y points and find
            % where the mask is defined
            z_line = reshape(cadat.STL.maskf(cySTLi,cxSTLi,:),[size(cadat.STL.maskf,3),1]);
            zs = find(isnan(z_line(1:czSTLi)),1,'last') + 1;
            ze = find(isnan(z_line(czSTLi:end)),1,'first') + czSTLi - 2;
            % Ensure zs and ze are defined
            if isempty(zs), zs = 1; end
            if isempty(ze), ze = length(x_line); end
            % Get the x line corresponding to the z and y points and find
            % where the mask is defined
            x_line = reshape(cadat.STL.maskf(cySTLi,:,czSTLi),[size(cadat.STL.maskf,2),1]);
            xs = find(isnan(x_line(1:cxSTLi)),1,'last') + 1;
            xe = find(isnan(x_line(cxSTLi:end)),1,'first') + cxSTLi - 2;
            % Ensure xs and xe are defined
            if isempty(xs), xs = 1; end
            if isempty(xe), xe = length(x_line); end
            % Get the y line corresponding to the x and z points and find
            % where the mask is defined
            y_line = reshape(cadat.STL.maskf(:,cxSTLi,czSTLi),[size(cadat.STL.maskf,1),1]);
            ys = find(isnan(y_line(1:cySTLi)),1,'last') + 1;
            ye = find(isnan(y_line(cySTLi:end)),1,'first') + cySTLi - 2;
            % Ensure ys and ye are defined
            if isempty(ys), ys = 1; end
            if isempty(ye), ye = length(y_line); end
            
            %%% USE VELOCITY MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get the z line corresponding to the x and y point and find
            % where the point is defined
% %             z_line = reshape(velmask(cyi,cxi,:),[size(velmask,3),1]);
% %             %zs = find(~isnan(z_line),1,'first');
% %             %ze = find(~isnan(z_line),1,'last');
% %             zs = find(isnan(z_line(1:czi)),1,'last') + 1;
% %             ze = find(isnan(z_line(czi:end)),1,'first') + czi - 2;
% %             %czi = round((ze+zs)/2);
% %             % Get the x and y lines corresponding to the x/y and z point.
% %             % Determine if the ellipse should start in the x-z or y-z plane
% %             x_line = reshape(velmask(cyi,:,czi),[size(velmask,2),1]);
% %             xs = find(isnan(x_line(1:cxi)),1,'last') + 1;
% %             xe = find(isnan(x_line(cxi:end)),1,'first') + cxi - 2;
% %             % Ensure xs and xe are defined
% %             if isempty(xs), xs = 1; end
% %             if isempty(xe), xe = length(x_line); end
% %             % Get y line and find end points
% %             y_line = reshape(velmask(:,cxi,czi),[size(velmask,1),1]);
% %             ys = find(isnan(y_line(1:cyi)),1,'last') + 1;
% %             ye = find(isnan(y_line(cyi:end)),1,'first') + cyi - 2;
% %             % Ensure ys and ye are defined
% %             if isempty(ys), ys = 1; end
% %             if isempty(ye), ye = length(y_line); end
            %%% USE VELOCITY MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Get x and y and z distances
            x_dist = xe - xs;
            y_dist = ye - ys;
            z_dist = ze - zs;
            
            % Want to use the smaller distance as the start for the
            % ellipse, then adjust it to be perpendicular to the center
            % line at that point
            if (x_dist <= y_dist) && (z_dist <= y_dist) % Want to use an x-z ellipse, with the center line moving primarily in the y-direction
                % Set the dx for the cross section
                dx1 = dx;
                dx2 = dz;
                % Get the ellipse point
                ell_pt = cyi;
                % Collect all points contained on the ellipse point plane
                ell_plane = permute(velmask(ell_pt,:,:),[2 3 1]);
                cc_points = ConnectedComponents(ell_plane,0,cxi,czi,ones(size(ell_plane)));
                % Set the diameter vector                
                diam_vec = [0,1,0];
                % Reshape indices to gather all in the ellipse point plane
                [zpl,xpl] = meshgrid(1:size(ell_plane,2),1:size(ell_plane,1));
                ypl = cyi*ones(size(xpl));
                indX = reshape(xpl,[size(xpl,1)*size(xpl,2),1]);
                indY = reshape(ypl,[size(ypl,1)*size(ypl,2),1]);
                indZ = reshape(zpl,[size(zpl,1)*size(zpl,2),1]);
                ccp = reshape(cc_points,[size(zpl,1)*size(zpl,2),1]);
                plX = indX(ccp > 0);
                plY = indY(ccp > 0);
                plZ = indZ(ccp > 0);             
                
                % Find the centerline to determine the normal of the vessel
                % Use a 7-plane stencil to minimize noise
                % This is done using the STL mask to improve accuracy
                % Find the ellipse plane in the STL mask closest to the
                % point of interest
                [~,cxSTLi] = min(abs(cadat.STL.x-cx));
                [~,cySTLi] = min(abs(cadat.STL.y-cy));
                [~,czSTLi] = min(abs(cadat.STL.z-cz));
                ct = 1; cntX = 0; cntY = 0; cntZ = 0; 
                for plct = -3:1:3
                    ell_cpl = permute(cadat.STL.maskf(cySTLi+plct,:,:),[2 3 1]);
                    cc_points = ConnectedComponents(ell_cpl,0,cxSTLi,czSTLi,ones(size(ell_cpl)));
                    BW = cc_points > 0;
                    cnt = regionprops(BW,'centroid');
                    cntX(ct,1) = cnt.Centroid(2);
                    cntZ(ct,1) = cnt.Centroid(1);
                    cntY(ct,1) = plct;
                    ct = ct + 1;
                    % Compute the diameter for the current STL plane
                    if plct == 0
                        bwstr = regionprops(BW,'area');
                        xs_area = bwstr.Area * cadat.STL.dx * cadat.STL.dz; % Cross sectional area in mm^2
                        bwstr = regionprops(BW,'EquivDiameter');
                        xs_diam = bwstr.EquivDiameter * cadat.STL.dx; % Diameter in mm
                    end
                end
                % Compute linear fit of center points to find angle of
                % normal
                ymag = 1;
                p = polyfit(cntY,cntX,1);
                xmag = p(1);
                p = polyfit(cntY,cntZ,1);
                zmag = p(1);
            elseif (y_dist <= x_dist) && (z_dist <= x_dist) % Want to use a y-z ellipse, with the center line moving primarily in the x-direction
                % Set the dx for the cross section
                dx1 = dz;
                dx2 = dy;
                % Get the ellipse point
                ell_pt = cxi;
                % Collect all points contained on the ellipse point plane
                ell_plane = permute(velmask(:,cxi,:),[1 3 2]);
                cc_points = ConnectedComponents(ell_plane,0,cyi,czi,ones(size(ell_plane)));
                % Set the diameter vector
                diam_vec = [1,0,0];
                % Reshape indices to gather all in the ellipse point plane
                [zpl,ypl] = meshgrid(1:size(ell_plane,2),1:size(ell_plane,1));
                xpl = cxi*ones(size(ypl));
                indX = reshape(xpl,[size(xpl,1)*size(xpl,2),1]);
                indY = reshape(ypl,[size(ypl,1)*size(ypl,2),1]);
                indZ = reshape(zpl,[size(zpl,1)*size(zpl,2),1]);
                ccp = reshape(cc_points,[size(zpl,1)*size(zpl,2),1]);
                plX = indX(ccp > 0);
                plY = indY(ccp > 0);
                plZ = indZ(ccp > 0);             
                
                % Find the centerline to determine the normal of the vessel
                % Use a 5-plane stencil to minimize noise
                [~,cxSTLi] = min(abs(cadat.STL.x-cx));
                [~,cySTLi] = min(abs(cadat.STL.y-cy));
                [~,czSTLi] = min(abs(cadat.STL.z-cz));
                ct = 1; cntX = 0; cntY = 0; cntZ = 0; 
                for plct = -3:1:3
                    ell_cpl = permute(cadat.STL.maskf(:,cxSTLi+plct,:),[1 3 2]);
                    cc_points = ConnectedComponents(ell_cpl,0,cySTLi,czSTLi,ones(size(ell_cpl)));
                    BW = cc_points > 0;
                    cnt = regionprops(BW,'centroid');
                    cntX(ct,1) = plct;
                    cntZ(ct,1) = cnt.Centroid(1);
                    cntY(ct,1) = cnt.Centroid(2);
                    ct = ct + 1;
                    % Compute the diameter for the current STL plane
                    if plct == 0
                        bwstr = regionprops(BW,'area');
                        xs_area = bwstr.Area * cadat.STL.dy * cadat.STL.dz; % Cross sectional area in mm^2
                        bwstr = regionprops(BW,'EquivDiameter');
                        xs_diam = bwstr.EquivDiameter * cadat.STL.dy; % Diameter in mm
                    end
                end
                % Compute linear fit of center points to find angle of
                % normal
                xmag = 1;
                p = polyfit(cntX,cntY,1);
                ymag = p(1);
                p = polyfit(cntX,cntZ,1);
                zmag = p(1);
            else % Want to use a x-y ellipse, with the center line moving primarily in the z-direction
                % Set the dx for the cross section
                dx1 = dx;
                dx2 = dy;
                % Get the ellipse point
                ell_pt = czi;
                % Collect all points contained on the ellipse point plane
                ell_plane = permute(velmask(:,:,czi),[1 2 3]);
                cc_points = ConnectedComponents(ell_plane,0,cyi,cxi,ones(size(ell_plane)));
                % Compute area and diameter of cross section
                BW = cc_points > 0;
                diam_vec = [0,0,1];
                % Reshape indices to gather all in the ellipse point plane
                [xpl,ypl] = meshgrid(1:size(ell_plane,2),1:size(ell_plane,1));
                zpl = czi*ones(size(ypl));
                indX = reshape(xpl,[size(xpl,1)*size(xpl,2),1]);
                indY = reshape(ypl,[size(ypl,1)*size(ypl,2),1]);
                indZ = reshape(zpl,[size(zpl,1)*size(zpl,2),1]);
                ccp = reshape(cc_points,[size(zpl,1)*size(zpl,2),1]);
                plX = indX(ccp > 0);
                plY = indY(ccp > 0);
                plZ = indZ(ccp > 0);             
                
                % Find the centerline to determine the normal of the vessel
                % Use a 5-plane stencil to minimize noise
                [~,cxSTLi] = min(abs(cadat.STL.x-cx));
                [~,cySTLi] = min(abs(cadat.STL.y-cy));
                [~,czSTLi] = min(abs(cadat.STL.z-cz));
                ct = 1; cntX = 0; cntY = 0; cntZ = 0; 
                for plct = -3:1:3
                    ell_cpl = permute(cadat.STL.maskf(:,:,czSTLi+plct),[1 2 3]);
                    cc_points = ConnectedComponents(ell_cpl,0,cySTLi,cxSTLi,ones(size(ell_cpl)));
                    BW = cc_points > 0;
                    cnt = regionprops(BW,'centroid');
                    cntZ(ct,1) = plct;
                    cntY(ct,1) = cnt.Centroid(1);
                    cntX(ct,1) = cnt.Centroid(2);
                    ct = ct + 1;
                    if plct == 0
                        bwstr = regionprops(BW,'area');
                        xs_area = bwstr.Area * cadat.STL.dy * cadat.STL.dx; % Cross sectional area in mm^2
                        bwstr = regionprops(BW,'EquivDiameter');
                        xs_diam = bwstr.EquivDiameter * cadat.STL.dz; % Diameter in mm
                    end
                end
                % Compute linear fit of center points to find angle of
                % normal
                zmag = 1;
                p = polyfit(cntZ,cntX,1);
                xmag = p(1);
                p = polyfit(cntZ,cntY,1);
                ymag = p(1);
            end
            % Make the normal a unit normal
            normal = [xmag,ymag,zmag];
            nmag = sqrt(sum(normal.^2));
            unormal = normal/nmag;
            unormal = diam_vec;
            % Adjust the diameter based on the angle of the unit normal
            costheta = sum(diam_vec.*unormal);
            diam_adj = xs_diam*costheta;
            area_adj = (pi/4)*diam_adj^2;
            
            % Save the cross section
            cadat.DATA.cross_section{q,1}.inds = [plX,plY,plZ];
            cadat.DATA.cross_section{q,1}.normal = unormal;
            cadat.DATA.cross_section{q,1}.diameter = diam_adj;
            cadat.DATA.cross_section{q,1}.area = area_adj;
            cadat.DATA.cross_section{q,1}.dx = [dx1,dx2];
        end
        cross_sections = cadat.DATA.cross_section;
        
        if cadat.options.SAVEOUTPUTS
            % Save cross sections
            save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'cross_sections.mat'],'cross_sections','xsecX','xsecY','xsecZ')
        end
    else
        % Load cross section
        fp = load([cadat.DIR.saveppfiles,cadat.FILE.anyname,'cross_sections.mat']);
        cadat.DATA.cross_section = fp.cross_sections;
        num_xsecs = size(cadat.DATA.cross_section,1);
    end
    
    % Initialize output
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    avg_vel = zeros(num_t,num_xsecs);
    re_num = zeros(num_t,num_xsecs);
    Q = zeros(num_t,num_xsecs);
    
    % Iterate through all time points
    fprintf('\nComputing flow rate...')
    prev_mark = 0;
    iter = 1;
    for q = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        % Load velocity file
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(q,'%05i'),'.mat']);
        % Get current velocities
        x = fp.x;
        y = fp.y;
        z = fp.z;
        u = fp.u;
        v = fp.v;
        w = fp.w;
        
        % Place the velocities on a grid
        cu = zeros(size(velmask));
        cv = zeros(size(velmask));
        cw = zeros(size(velmask));
        for ptn = 1:1:length(x)
            [~,cind] = min(abs(maskX - x(ptn)));
            [~,rind] = min(abs(maskY - y(ptn)));
            [~,zind] = min(abs(maskZ - z(ptn)));
            
            % Put velocity in corresponding spatial location
            cu(rind,cind,zind) = u(ptn);
            cv(rind,cind,zind) = v(ptn);
            cw(rind,cind,zind) = w(ptn);
        end
        
        % Iterate through each cross section to identify the flow rate at
        % that time
        for xc = 1:1:num_xsecs
            % Get current cross section data
            currpts = cadat.DATA.cross_section{xc,1}.inds;
            currnorm = cadat.DATA.cross_section{xc,1}.normal;
            %currnorm = currnorm/(sqrt(sum(currnorm.^2))); % Make unit normal
            currarea = cadat.DATA.cross_section{xc,1}.area; % In mm^2
            currdiam = cadat.DATA.cross_section{xc,1}.diameter; % In mm
            currdx = cadat.DATA.cross_section{xc,1}.dx; % In mm
            num_pts = size(currpts,1);
            %currdx(1) = sqrt(currarea/num_pts);
            %currdx(2) = currdx(1);
            % Iterate through all points in the cross section
            Qf = 0; % Initialize final flow
            Qf_2c = 0; % Initialize 2-component flow rate
            for pts = 1:1:size(currpts,1)
                % Get velocities at current point (velocities in m/s)
                curr_u = cu(currpts(pts,2),currpts(pts,1),currpts(pts,3));
                curr_v = cv(currpts(pts,2),currpts(pts,1),currpts(pts,3));
                curr_w = cw(currpts(pts,2),currpts(pts,1),currpts(pts,3));
                
                % Multiply velocities by normal
                norm_u = curr_u*currnorm(1);
                norm_v = curr_v*currnorm(2);
                norm_w = curr_w*currnorm(3);
                
                % Compute flow rate addition, scale appropriately by
                % integral area
                Qc = (norm_u + norm_v + norm_w)*(currdx(1)*1e-3*currdx(2)*1e-3);
                Qf = Qf + Qc; % Final flow rate
                Qsave(pts,1) = Qc;
                save_vel(pts,:) = [curr_u,curr_v,curr_w];
                %Qc_2c = (norm_u + norm_v)*(dx*1e-3*dz*1e-3);
                %Qf_2c = Qf_2c + Qc_2c;
            end
            Qf = Qf;%*currarea*(1e-3)^2/num_pts;
            % Save data for current cross section at current time step
            Q(iter,xc) = Qf; % In m^3/s
            %Q_2c(iter,xc) = Qf_2c; % In m^3/s
            avg_vel(iter,xc) = Qf/(currarea*(1e-3)^2); % Average velocity computed from flow rate (in m/s)
            re_num(iter,xc) = (currdiam*1e-3)*abs(avg_vel(iter,xc))/cadat.DATA.nu; % Current reynolds number
        end
        % Print the percent completed for steps of 10
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
        
        % Iterate counter
        iter = iter + 1;
    end
    % Store final data
    cadat.DATA.flow_rate.Q = Q; % In m^3/s
    %cadat.DATA.flow_rate.Q_2c = Q_2c;
    cadat.DATA.flow_rate.Re = re_num;
    cadat.DATA.flow_rate.avgvel = avg_vel; % In m/s
    
    if cadat.options.SAVEOUTPUTS
        % Save flow rate data
        flow_rate = cadat.DATA.flow_rate;
        save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'flow_rate_data.mat'],'flow_rate')
        fprintf(' completed successfully\n')
    end
%else
    % Load flow rate data
%    fp = load([cadat.DIR.saveppfiles,cadat.FILE.anyname,'flow_rate_data.mat']);
%    cadat.DATA.flow_rate = fp.flow_rate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SECONDARY FLOW VARIABLE CALCULATION OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE WALL NORMAL VECTORS
if cadat.options.COMPNORMS
    fprintf('\n----- Wall Normal Calculation -----')
    % Load the velocity mask
    %fp = load([cadat.DIR.savefiles,cadat.FILE.maskgrid]);
    %velmask = fp.velmask;
    xmask = maskX; %fp.xmask;
    ymask = maskY; %fp.ymask;
    zmask = maskZ; %fp.zmask;

    % Get the STL mask
    stlmask = cadat.STL.maskf;
    [Xstl,Ystl,Zstl] = meshgrid(cadat.STL.x,cadat.STL.y,cadat.STL.z);
    % Set parameters
    smooth_norms = cadat.options.SMOOTHNORMS;
    plot_norms = 0;%cadat.options.DEBUGPLOTS;
    
    % Compute edge normals
    if strcmp(cadat.FILE.filetype,'MRI') || cadat.options.USEVOXAVG == 1
        % Load the appropriate edge file based on the modality being used
        if strcmp(cadat.FILE.filetype,'MRI') || strcmp(cadat.FILE.filetype,'CFD')
            % Load the edge normals
            loadname = ['any',cadat.FILE.anynum,'_','CFD','_full_resolution_'];
            fp = load([cadat.DIR.saveppfiles,loadname,'edge_normals_smoothed.mat']);
            highres_norm = fp.edge_normal;
            % Load the X,Y,Z locations of the high resolution normals
            loadgrid = ['ANY-',cadat.FILE.anynum,'_','CFD','_masked_velocity_grid'];
            fp = load([cadat.DIR.savefiles,loadgrid,'.mat']);
            x_highres = fp.xmask;
            y_highres = fp.ymask;
            z_highres = fp.zmask;
        else
            % Load the edge normals
            loadname = ['any',cadat.FILE.anynum,'_','STB','_full_resolution_'];
            fp = load([cadat.DIR.saveppfiles,loadname,'edge_normals_smoothed.mat']);
            highres_norm = fp.edge_normal;
            % Load the X,Y,Z locations of the high resolution normals
            loadgrid = ['ANY-',cadat.FILE.anynum,'_','STB','_masked_velocity_grid'];
            fp = load([cadat.DIR.savefiles,loadgrid,'.mat']);
            x_highres = fp.xmask;
            y_highres = fp.ymask;
            z_highres = fp.zmask;
        end
        % Form grid of high resolution data
        [Xhr,Yhr,Zhr] = meshgrid(x_highres,y_highres,z_highres);
        % Reshape grid to lines
        Xhr_l = reshape(Xhr,[size(Xhr,1)*size(Xhr,2)*size(Xhr,3),1]);
        Yhr_l = reshape(Yhr,[size(Xhr,1)*size(Xhr,2)*size(Xhr,3),1]);
        Zhr_l = reshape(Zhr,[size(Xhr,1)*size(Xhr,2)*size(Xhr,3),1]);
        highres_norm_l = reshape(highres_norm,[size(Xhr,1)*size(Xhr,2)*size(Xhr,3),1]);
        % Keep only the points whose edges are non-zero
        xhr_edge = []; yhr_edge = []; zhr_edge = []; clear hr_edge
        ct = 1;
        for zz = 1:1:size(Xhr_l)
            cnorm = highres_norm_l{zz,1};
            if ~isempty(cnorm)
                xhr_edge(ct,1) = Xhr_l(zz);
                yhr_edge(ct,1) = Yhr_l(zz);
                zhr_edge(ct,1) = Zhr_l(zz);
                hr_edge{ct,1} = cnorm;
                ct = ct + 1;
            end
        end
        
        % Find the edges in the current velocity field
        [Xm,Ym,Zm] = meshgrid(xmask,ymask,zmask);
        edge_pts = zeros(size(Zm));
        for ii = 1:1:size(Xm,1)
            for jj = 1:1:size(Xm,2)
                for kk = 1:1:size(Xm,3)
                    % Determine if the current point is in the boundary
                    if ~isnan(velmask(ii,jj,kk))
                        % Determine if point is an edge
                        % Check if point is at velocity mask edge
                        ii_rng = (ii-1 > 0) && (ii+1 < size(Xm,1));
                        jj_rng = (jj-1 > 0) && (jj+1 < size(Xm,2));
                        kk_rng = (kk-1 > 0) && (kk+1 < size(Xm,3));
                        if ii_rng && jj_rng && kk_rng
                            x_neighs = ~isnan(velmask(ii-1,jj,kk)) + ~isnan(velmask(ii+1,jj,kk));     
                            y_neighs = ~isnan(velmask(ii,jj-1,kk)) + ~isnan(velmask(ii,jj+1,kk));
                            z_neighs = ~isnan(velmask(ii-1,jj,kk-1)) + ~isnan(velmask(ii,jj,kk+1));
                            all_neighs_gd = x_neighs + y_neighs + z_neighs;
                            if all_neighs_gd < 6
                                % Set the edge point to 1
                                edge_pts(ii,jj,kk) = 1;
                            end       
                        else
                            edge_pts(ii,jj,kk) = 1;
                        end
                    end
                end
            end
        end
        % Obtain all boundary points as a line -  for plotting purposes
        x_all_pts = reshape(Xm,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
        y_all_pts = reshape(Ym,[size(Ym,1)*size(Ym,2)*size(Ym,3),1]);
        z_all_pts = reshape(Zm,[size(Zm,1)*size(Zm,2)*size(Zm,3),1]);
        is_edge = reshape(edge_pts,[size(Xm,1)*size(Xm,2)*size(Xm,3),1]);
        % Limit to the points in the mask only
        x_epts = x_all_pts(is_edge == 1);%(~isnan(mask_pts));
        y_epts = y_all_pts(is_edge == 1); %~isnan(mask_pts));
        z_epts = z_all_pts(is_edge == 1); %~isnan(mask_pts));
        
        % For each edge normal obtain the edge on the full resolution map
        % that is closest to the edge point
        edge_normal = cell(size(edge_pts)); % Initialize edge points
        for ii = 1:1:size(Xm,1)
            for jj = 1:1:size(Xm,2)
                for kk = 1:1:size(Xm,3)
                    % Only calculate if the point is an edge
                    if edge_pts(ii,jj,kk) == 1
                        %%% Get current edge point info %%%%%%%%%%%%%%%%%%%
                        x = Xm(ii,jj,kk);
                        y = Ym(ii,jj,kk);
                        z = Zm(ii,jj,kk);
                        
                        % Find the closest edge point on the high
                        % resolution edge data
                        pt_dist = sqrt((xhr_edge-x).^2+(yhr_edge-y).^2+(zhr_edge-z).^2);
                        [min_dist,min_ind] = min(pt_dist);
                        opt_norm = hr_edge{min_ind,1};
                        % Set the edge normal to the closest edge
                        edge_normal{ii,jj,kk} = opt_norm;
                    end
                end
            end
        end
        % Plot the normals
        % Get all of the normals into arrays
        Unorm = zeros(size(edge_normal));
        Vnorm = zeros(size(edge_normal));
        Wnorm = zeros(size(edge_normal));
        % Iterate through all points and separate the normals
        for ii = 1:1:size(Unorm,1)
            for jj = 1:1:size(Unorm,2)
                for kk = 1:1:size(Unorm,3)
                    curr_norm = edge_normal{ii,jj,kk};
                    if ~isempty(curr_norm)
                        Unorm(ii,jj,kk) = curr_norm(1);
                        Vnorm(ii,jj,kk) = curr_norm(2);
                        Wnorm(ii,jj,kk) = curr_norm(3);
                    end
                end
            end
        end
        % Plot the normals
        figure(27); hold off; scatter3(x_epts,y_epts,z_epts,3,'MarkerEdgeColor',[0 0.447 0.741],'MarkerFaceColor',[0 0.447 0.741]); hold on;
        quiver3(Xm,Ym,Zm,Unorm,Vnorm,Wnorm,'r');
        
    else
        edge_normal = compute_edge_normals_v2(velmask,xmask,ymask,zmask,stlmask,Xstl,Ystl,Zstl,smooth_norms,plot_norms);
    end
    
    if cadat.options.SAVEOUTPUTS
        % Save the edge normals
        save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'edge_normals_smoothed.mat'],'edge_normal');
    end
    fprintf('-----------------------------------\n')
elseif cadat.options.LOADNORMS
    % Load the smoothed edge normals
    fp = load([cadat.DIR.saveppfiles,cadat.FILE.anyname,'edge_normals_smoothed.mat']);
    edge_normal = fp.edge_normal;
end


%% COMPUTE VELOCITY GRADIENTS
% Gradients are computed using 3D radial basis functions - thin plate
% splines used for interpolation of the RBF
if cadat.options.COMPVELGRAD
    addpath /Users/Melissa/Documents/Code/general_codes/gradients
    fprintf('\n----- Computing gradients ----- \n')
    % Set window size for wall shear stress
    winsize = 3;
    
    % Convert the grid to meters (initially it is in mm)
    Vxm = Vx*1e-3;
    Vym = Vy*1e-3;
    Vzm = Vz*1e-3;
    
    % Iterate through all time points
    for curr_time = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        fprintf('File %i: ',curr_time)
        % Load the velocity of the current time step
        fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(curr_time,'%05i'),'.mat']);
        % Get current velocities
        x = fp.x;
        y = fp.y;
        z = fp.z;
        u = fp.u;
        v = fp.v;
        w = fp.w;
        % Place the velocities on a grid
        Vu = zeros(size(velmask));
        Vv = zeros(size(velmask));
        Vw = zeros(size(velmask));
        for ptn = 1:1:length(x)
            [~,cind] = min(abs(maskX - x(ptn)));
            [~,rind] = min(abs(maskY - y(ptn)));
            [~,zind] = min(abs(maskZ - z(ptn)));  
            % Put velocity in corresponding spatial location
            Vu(rind,cind,zind) = u(ptn);
            Vv(rind,cind,zind) = v(ptn);
            Vw(rind,cind,zind) = w(ptn);
        end
        
        % Compute the spatial gradients for the current time step
        [dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz] = gradient_rbf_3D(Vxm,Vym,Vzm,Vu,Vv,Vw,velmask,winsize);
        
        if cadat.options.SAVEOUTPUTS
            % Save the edge normals
            save([cadat.DIR.saveppfiles,'velocity_gradients',fslash,cadat.FILE.anyname,'velgrads_',num2str(curr_time,'%05i'),'.mat'],'dudx','dudy','dudz','dvdx','dvdy','dvdz','dwdx','dwdy','dwdz');
        end
    end
    fprintf(' complete\n')
end

%% COMPUTE DIVERGENCE
% Computes divergence using the velocity gradients previously computed
if cadat.options.COMPDIVERGENCE
    fprintf('\nComputing divergence...')
    % Iterate through all time
    iter = 1; prev_mark = 0; 
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    avg_div_3c = zeros(num_t,1);
    avg_div_2c = zeros(num_t,1);
    med_div_3c = zeros(num_t,1);
    med_div_2c = zeros(num_t,1);
    for curr_time = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        % Load the gradients for the current file
        fp = load([cadat.DIR.saveppfiles,'velocity_gradients',fslash,cadat.FILE.anyname,'velgrads_',num2str(curr_time,'%05i'),'.mat']);
        % Compute the divergence
        div_3c = fp.dudx + fp.dvdy + fp.dwdz;
        div_3c(1:10,1:10,:) = 0;
        div_2c = fp.dudx + fp.dvdy;
        div_2c(1:10,1:10,:) = 0;
        div_3c_line = abs(div_3c(:));
        div_3c_line = div_3c_line(div_3c_line > 0);
        div_2c_line = abs(div_2c(:));
        div_2c_line = div_2c_line(div_2c_line > 0);
        avg_div_3c(iter,1) = mean(div_3c_line);
        avg_div_2c(iter,1) = mean(div_2c_line);
        med_div_3c(iter,1) = median(div_3c_line);
        med_div_2c(iter,1) = median(div_2c_line);
        if cadat.options.SAVEOUTPUTS
            % Save the divergence field
            save([cadat.DIR.saveppfiles,'divergence_fields',fslash,cadat.FILE.anyname,'divergence_',num2str(curr_time,'%05i'),'.mat'],'div_3c','div_2c');
        end
        
        % Print the percent completed for steps of 10
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
        
        % Update iteration counter
        iter = iter + 1;
    end
    % Save the average divergence time series
    save([cadat.DIR.saveppfiles,'divergence_fields',fslash,cadat.FILE.anyname,'divergence_time_series.mat'])
    fprintf(' complete\n')
end


%% COMPUTE WALL SHEAR STRESS
% Computes the wall shear stress using the velocity graidents and wall
% normals previously computed
if cadat.options.COMPWSS
    grid_size = 3;
    fprintf('\nComputing wall shear stress...')
    % Iterate through all time
    iter = 1; prev_mark = 0; 
    num_t = length(cadat.options.TS:cadat.options.Tskip:cadat.options.TE);
    for curr_time = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        % Load the velocity gradients for the current time step
        fp = load([cadat.DIR.saveppfiles,'velocity_gradients',fslash,cadat.FILE.anyname,'velgrads_',num2str(curr_time,'%05i'),'.mat']);
        dudx = fp.dudx; dudy = fp.dudy; dudz = fp.dudz;
        dvdx = fp.dvdx; dvdy = fp.dvdy; dvdz = fp.dvdz;
        dwdx = fp.dwdx; dwdy = fp.dwdy; dwdz = fp.dwdz;
        
        % Compute wall shear stress (in Pa)
        [wss_mag,wss_x,wss_y,wss_z] = compute_wss(edge_normal,cadat.DATA.mu,grid_size,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz);
        
        % Save the wall shear stress output
        if cadat.options.SAVEOUTPUTS
            save([cadat.DIR.saveppfiles,'wall_shear_stress',fslash,cadat.FILE.anyname,'wss_',num2str(curr_time,'%05i'),'.mat'],'wss_mag','wss_x','wss_y','wss_z')
        end
        
        % Print the percent completed for steps of 10
        perc_complete = floor(100*iter/num_t);
        if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
            fprintf(' %i%%',perc_complete)
            prev_mark = perc_complete;
        end
        
        % Increment the iteration number
        iter = iter + 1;
    end
    
    % Compute the time averaged wall shear stress
    tawss_x = zeros(size(wss_x));
    tawss_y = zeros(size(wss_y));
    tawss_z = zeros(size(wss_z));
    tawss_mag = zeros(size(wss_mag));
    for curr_time = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        % Load the current wall shear stress 
        fp = load([cadat.DIR.saveppfiles,'wall_shear_stress',fslash,cadat.FILE.anyname,'wss_',num2str(curr_time,'%05i'),'.mat']);
        % Add the wall shear stress magnitude to the time averaged array
        tawss_x = tawss_x + abs(fp.wss_x);
        tawss_y = tawss_y + abs(fp.wss_y);
        tawss_z = tawss_z + abs(fp.wss_z);
        tawss_mag = tawss_mag + fp.wss_mag;
    end
    % Average the wall shear stress
    tawss_x = tawss_x/num_t;
    tawss_y = tawss_y/num_t;
    tawss_z = tawss_z/num_t;
    tawss_mag = tawss_mag/num_t;
    if cadat.options.SAVEOUTPUTS
        save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'TAWSS','.mat'],'tawss_mag','tawss_x','tawss_y','tawss_z')
    end
    fprintf(' complete\n')
end

%% COMPUTE OSCILLATORY SHEAR INDEX (OSI)
if cadat.options.COMPOSI
    % Load the wall shear stress at each time step and sum according to the
    % OSI formula - this is done for x, y, and z wall shear stress
    for curr_time = cadat.options.TS:cadat.options.Tskip:cadat.options.TE
        fp = load([cadat.DIR.saveppfiles,'wall_shear_stress',fslash,cadat.FILE.anyname,'wss_',num2str(curr_time,'%05i'),'.mat']);
        if curr_time == cadat.options.TS
            % For first time step, initialize the average vectors
            wss_x_avg = zeros(size(fp.wss_x));
            wss_y_avg = zeros(size(fp.wss_y));
            wss_z_avg = zeros(size(fp.wss_z));  
            wss_x_absavg = zeros(size(fp.wss_x));
            wss_y_absavg = zeros(size(fp.wss_y));
            wss_z_absavg = zeros(size(fp.wss_z));   
        end
        % Add the current values to the average vectors
        wss_x_avg = wss_x_avg + fp.wss_x;
        wss_y_avg = wss_y_avg + fp.wss_y;
        wss_z_avg = wss_z_avg + fp.wss_z;
        wss_x_absavg = wss_x_absavg + abs(fp.wss_x);
        wss_y_absavg = wss_y_absavg + abs(fp.wss_y);
        wss_z_absavg = wss_z_absavg + abs(fp.wss_z);
    end
    
    % Compute the OSI for each coordinate
    osi_x = 0.5*(1 - abs(wss_x_avg)./wss_x_absavg);
    osi_y = 0.5*(1 - abs(wss_y_avg)./wss_y_absavg);
    osi_z = 0.5*(1 - abs(wss_z_avg)./wss_z_absavg);
    osi_mag = sqrt(osi_x.^2 + osi_y.^2 + osi_z.^2);
    osi_avg = (osi_x + osi_y + osi_z)/3;
    % Save the wall shear stress output
    if cadat.options.SAVEOUTPUTS
        save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'OSI','.mat'],'osi_mag','osi_x','osi_y','osi_z','osi_avg')
    end
    
end

%% COMPUTE RELATIVE RESIDENCE TIMES (RRT)
% RRT is computed as: RRT = 1/((1-2*OSI)*TAWSS)
if cadat.options.COMPRRT
    % Load TAWSS
    fp = load([cadat.DIR.saveppfiles,cadat.FILE.anyname,'TAWSS','.mat']);
    tawss_x = fp.tawss_x;
    tawss_y = fp.tawss_y;
    tawss_z = fp.tawss_z;
    tawss_mag = fp.tawss_mag;
    % Load OSI
    fp = load([cadat.DIR.saveppfiles,cadat.FILE.anyname,'OSI','.mat']);
    osi_avg = fp.osi_avg; % Use coordinate averaged osi
    osi_x = fp.osi_x;
    osi_y = fp.osi_y;
    osi_z = fp.osi_z;
    
    % Compute RRT
    rrt_x = 1./((1-2*osi_x).*tawss_x);
    rrt_y = 1./((1-2*osi_y).*tawss_y);
    rrt_z = 1./((1-2*osi_z).*tawss_z);
    rrt_avg = (rrt_x + rrt_y + rrt_z)/3;
    
    % Save the wall shear stress output
    if cadat.options.SAVEOUTPUTS
        save([cadat.DIR.saveppfiles,cadat.FILE.anyname,'RRT','.mat'],'rrt_avg','rrt_x','rrt_y','rrt_z')
    end
    
end

%% COMPUTE VORTICITY
% Computes the vorticity using the velocity graidents previously computed
if cadat.options.COMPVORT
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTTING OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VELOCITY STREAMLINE PLOTS
if cadat.options.PLOTSTREAMTRACES
    plot_vel_streamtraces(cadat);
end

%% VELOCITY POINT COMPARISON PLOTS
% This step uses the velocity fields and cross sections computed
if cadat.options.PLOTVELCOMP
    plot_vel_comp(cadat);
end


%% VELOCITY PROFILE PLOTS
% This step uses the velocity fields and cross sections computed
if cadat.options.PLOTVELPROFILES
    plot_vel_profiles(cadat,fslash,0);
end

%% DIVERGENCE PLOTS
if cadat.options.PLOTDIVERGENCE
    plot_divergence(cadat,fslash);
end

%% WALL SHEAR STRESS PLOTS
if cadat.options.PLOTWSS
    plot_wall_shear_stress(cadat,fslash);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APPENDIX - OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % Determine the points that are separated by NaN
                        % values (cross through invalid mask regions). Do
                        % this only for the first 500 points, which should
                        % provide a sufficient cohort of points
%                         neigh_valid = zeros(500,1);
%                         neigh_valid(1) = 1;
%                         for qq = 2:1:500
%                             % Get current neighbor information
%                             neighX = x_epts(sinds(qq));
%                             neighY = y_epts(sinds(qq));
%                             neighZ = z_epts(sinds(qq));
%                             % Get line of points from current point to
%                             % possible neighbor point
%                             xl = x + (neighX-x)*neigh_t;
%                             yl = y + (neighY-y)*neigh_t;
%                             zl = z + (neighZ-z)*neigh_t;
%                             % Limit the line to only between the two points
%                             minX = min(neighX,x); maxX = max(neighX,x);
%                             minY = min(neighY,y); maxY = max(neighY,y);
%                             minZ = min(neighZ,z); maxZ = max(neighZ,z);
%                             kpl = (xl>=minX) & (xl<=maxX) & (yl>=minY) & (yl<=maxY) & (zl>=minZ) & (xl<=maxZ);
%                             xl = xl(kpl);
%                             yl = yl(kpl);
%                             zl = zl(kpl);
%                             % Convert the line to pixels
%                             line_inds = [];
%                             for ptnum = 1:1:length(xl)
%                                 % Current line point
%                                 clx = xl(ptnum);
%                                 cly = yl(ptnum);
%                                 clz = zl(ptnum);
%                                 % Nearest mask point
%                                 if (clx>=min(xmask)) && (clx<=max(xmask)) && (cly>=min(ymask)) && (cly<=max(ymask)) && (clz>=min(zmask)) && (clz<=max(zmask))
%                                     [~,mx_ind] = min(abs(xmask-clx));
%                                     [~,my_ind] = min(abs(ymask-cly));
%                                     [~,mz_ind] = min(abs(zmask-clz));
%                                     line_inds = [line_inds;my_ind,mx_ind,mz_ind]; % Order as (row,col,z)
%                                 end
%                             end
%                             % Keep only the unique line index rows
%                             line_inds = unique(line_inds,'rows');
%                             line_pix = zeros(size(line_inds,1),1);
%                             for rownum = 1:1:size(line_inds,1)
%                                 line_pix(rownum) = velmask(line_inds(rownum,1),line_inds(rownum,2),line_inds(rownum,3));
%                             end
%                             % Determine if any invalid points exist in the
%                             % line
%                             if sum(isnan(line_pix)) == 0
%                                 neigh_valid(qq) = 1;
%                             end
%                             
%                         end
%                         if sum(neigh_valid(1:25)==0)>0
%                             keyboard
%                         end
%                         sinds_kp = sinds(neigh_valid==1);
%                         sinds_rej = sinds(neigh_valid==0);
%                         kp_inds = sinds_kp(1:25);
                        %ptsurfX = x_epts(kp_inds);
                        %ptsurfY = y_epts(kp_inds);
                        %ptsurfZ = z_epts(kp_inds);





