clear cadat

%% DIRECTORY AND FILE NAME INFORMATION
cadat.FILE.anynum = '007'; % Aneurysm identifiable number
cadat.FILE.filetype = 'STB'; % Should be STB (Shake the box), TOMO, MRI, or CFD
    % STB loads mat file, TOMO loads dat file
    % MRI loads dat file
cadat.SYSTEM.cpu = 'pc';
cadat.FILE.run_num = '1';

% Set directory info based on system
if strcmp(cadat.SYSTEM.cpu,'pc')
    cadat.SYSTEM.file_base = 'Z:\Projects\Cerebral_Aneurysm\ANY\';
    fslash = '\';
else
    cadat.SYSTEM.file_base = '~/Projects/Cerebral_Aneurysm/ANY/';
    fslash = '/';
end

cadat.DIR.velfiles = [cadat.SYSTEM.file_base,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'raw',fslash];
cadat.DIR.velout = [cadat.SYSTEM.file_base,cadat.FILE.anynum,fslash,cadat.FILE.filetype,fslash,'vel_registered_masked',fslash];
cadat.DIR.stlfiles = [cadat.SYSTEM.file_base,cadat.FILE.anynum,fslash,'geometry',fslash];
cadat.DIR.savefiles = [cadat.SYSTEM.file_base,cadat.FILE.anynum,fslash,'registration_files',fslash];
if strcmp(cadat.FILE.filetype,'STB')
    cadat.FILE.basename = ['any',cadat.FILE.anynum,'_stb_v',cadat.FILE.run_num,'_'];
    cadat.FILE.baseinname = ['any',cadat.FILE.anynum,'_stb_vi_'];%,cadat.FILE.run_num,'_'];
elseif strcmp(cadat.FILE.filetype,'CFD')
    cadat.FILE.basename = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_'];
    cadat.FILE.baseinname = ['any',cadat.FILE.anynum,'_cfd_'];
else
    cadat.FILE.basename = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_'];
    cadat.FILE.baseinname = ['any',cadat.FILE.anynum,'_',cadat.FILE.filetype,'_'];
end

cadat.FILE.ext = '.mat'; % Not currently used

%% SPATIAL RESOLUTION
% The spatial resolution controls the uncertainty in the registration. A
% finer resolution means the registration is done with more precision (with
% a finer resolution). However the resolution needs to be big enough that
% there are enough points in the masks. The STL typically has many points
% so this is limited by the resolution of the data. For STB and CFD this
% can be pushed higher. For gridded data such as MRI and TOMO, the
% resolution likely needs to be larger.
% For TOMO, this likely needs to be the dx dy dz from those fields.
cadat.DATA.dx = 0.4; % in mm
cadat.DATA.dy = 0.4; % in mm
cadat.DATA.dz = 0.4; % in mm

% The STL typically has many points so it can use a higher resolution.
% However, the higher the STL resolution, the longer the STL mask 
% generation will be.
cadat.STL.dx = 0.2; % in mm
cadat.STL.dy = 0.2; % in mm
cadat.STL.dz = 0.2; % in mm

%% TIME STEP INFORMATION
cadat.TIME.ts = 1;
cadat.TIME.te = 5766;

%% VELOCITY SCALE AND FLIP INFORMATION
% Velocity files need to be in mm for x, y, z coordinates and m/s for u, v,
% and w velocities
% Set if the file is in mm or m
cadat.DATA.inMM = 1; % 1 - File is in mm, 0 - File is in m

%% STL INFORMATION
% STL file either needs to be imported and converted to a grid format, or
% the adjusted STL mat file needs to be loaded in. This step only needs to
% be completed one time so if multiple time points are to be analyzed, it
% will only run it for the 1st time point.
cadat.OPTIONS.ADJSTL = 0; % If STL file needs to be imported 
cadat.STL.fp_stl_name = 'NW-ANY-007_GEOM_SMOOTH_05MAY2016_FINAL_EXTADDED'; % DO NOT include extension
cadat.STL.stl_save_name = ['STL-ANY-',cadat.FILE.anynum];
cadat.STL.adjstlptsfile = [cadat.STL.stl_save_name,'_gridded_filtered'];
cadat.STL.adjstlflipfile = [cadat.STL.stl_save_name,'_gridded_filtered_flipped'];

%% STL AND DATA FLIP INFORMATION
cadat.OPTIONS.DETERMINE_FLIPS = 0; % Determine if the STL file or data, velocity fields need to be flipped 
    % This step requires user inputs
cadat.OPTIONS.save_stl_flips = 0; % Choose whether or not to save the STL file data since this only needs
    % to be done a single time across all datasets
cadat.DATA.flipfilename = [cadat.STL.stl_save_name,'_',cadat.FILE.filetype,'_flip-details.mat'];

%% STL MASK INFORMATION
% Creates mask for each z-plane from the adjusted STL file, or loads in the
% mask if it has previously been created. This step only needs to
% be completed one time so if multiple time points are to be analyzed, it
% will only run it for the 1st time point.
cadat.OPTIONS.CREATEMASK = 0; % Choose whether to create the STL mask
cadat.STL.maskfilename = [cadat.STL.stl_save_name,'_STLMask.mat'];

%% FILTERED/DATA MASK INFORMATION
% Creates mask of where non-zero velocity is, or loads in the mask if it
% has previously been created. This step only needs to be completed one
% time so if multiple time points are to be analyzed, it will only run it
% for the 1st time point.
cadat.OPTIONS.CREATEDATAMASK = 1; % Choose whether to create the data velocity mask or load it
cadat.VELMASK.velmaskfilename = [cadat.STL.stl_save_name,'_',cadat.FILE.filetype,'_DataMask.mat'];
cadat.VELMASK.ts = 1; % Starting time for creating the velocity mask
cadat.VELMASK.te = 10; % Ending time for creating the velocity mask
cadat.VELMASK.tstep = 1; % Time step for the velocity mask
cadat.STL.filtmaskname = [cadat.STL.stl_save_name,'_',cadat.FILE.filetype,'_STLMask_filt.mat'];
cadat.VELMASK.filtmaskname = [cadat.STL.stl_save_name,'_',cadat.FILE.filetype,'_DataMask_filt.mat'];
cadat.OPTIONS.TESTRESOLUTION = 0; % Choose wheter to just test the resolution for mask
    % If this step is on, all other options will not be processed and only
    % the data mask will be evaluated and the rest of the program will not
    % run.
    
%% REGISTRATION SHIFT INFORMATION
% Register the velocity fields to the STL by adjusting indices. This step
% only needs to be completed one time so if multiple time points are to be 
% analyzed, it will only run it for the 1st time point.
cadat.OPTIONS.FINDOPTSHIFT = 1; % Choose whether to find the optimal shift or load it
cadat.SHIFT.optshiftfile = [cadat.STL.stl_save_name,'_',cadat.FILE.filetype,'_optshift.mat'];
cadat.SHIFT.optrotfile = [cadat.STL.stl_save_name,'_',cadat.FILE.filetype,'_optrot.mat'];

%% VELOCITY REGISTRATION INFOMRATION
% Masks and registers the active velocity field (current time step)
cadat.OPTIONS.MASKREGVELOCITY = 1; % Choose whether to register and mask the velocity field

%% SAVE REGISTERED AND MASKED VELOCITY FIELD
% Mask and registering velocity **MUST** be turned **ON** for these to save
% correctly. Saving the zeros is ONLY recommended if TOMO is being used
% because it is on a grid. For all other methods (MRI, CFD, STB), do not
% save zeros because the data still needs to be converted to a grid
cadat.FILE.outputfilebase = [cadat.FILE.basename,'regmask_'];
cadat.OPTIONS.SAVE_RM_DATS = 1; % Choose whether to save the registered and masked velocity field
cadat.OPTIONS.SAVE_ZEROS_DAT = 0; % Choose whether to save the zero values for a uniform grid (DAT file only)
cadat.OPTIONS.SAVE_MAT = 1; % Save .mat registered and masked velocity files
cadat.OPTIONS.SAVE_ZEROS_MAT = 0; % Choose whether to save the zero values for a uniform grid (MAT file only)


%% RUN THE CODE
[cadat] = ca_register_vel_to_stl_v3(cadat);


