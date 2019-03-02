
dirpath = '~/ANY/007/MRI/';
outpath = '~/ANY/007/MRI/raw_unflipped/';
fmr = load([dirpath,'Vasculature_mask_struct.mat']);
vel = load([dirpath,'vel_struct.mat']);
out_file = 'any007_mri_';
ts = 1;
te = 13;

% Get the mask
mask = fmr.mrStruct.dataAy;
mask = permute(mask,[2 1 3]); % Permute to make y row and x col
mask(:,1:55,:) = 0;
mask(:,95:end,:) = 0;
mask(1:70,:,:) = 0;
mask(100:end,:,:) = 0;
mask(90:101,61:73,:) = 0;

% Get the velocity files
[n,m,k,c,t] = size(vel.mrStruct.dataAy);
mag = vel.mrStruct.vox;
xl = mag(1)*(1:1:n);
yl = mag(2)*(1:1:m);
zl = mag(3)*(1:1:k);
[X,Y,Z] = meshgrid(xl,yl,zl);

Vx = zeros(m,n,k,t);
Vy = zeros(m,n,k,t);
Vz = zeros(m,n,k,t);
for q = 1:1:t
    % Get current velocities
    cu = double(vel.mrStruct.dataAy(:,:,:,1,q));
    cv = double(vel.mrStruct.dataAy(:,:,:,2,q));
    cw = double(vel.mrStruct.dataAy(:,:,:,3,q));
    
    % Permute so row is y, col is x
    cu = permute(cu,[2 1 3]);
    cv = permute(cv,[2 1 3]);
    cw = permute(cw,[2 1 3]);
    
    % Save each time step
    Vx(:,:,:,q) = cu;
    Vy(:,:,:,q) = cv;
    Vz(:,:,:,q) = cw;
end

%%% Find trimmed indices from the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Easier to set by user!!! mask arrays computed for user to set trimmed
% indices.
mask_zsum = sum(mask,3);
mask_r = sum(mask_zsum,2);
mask_c = transpose(sum(mask_zsum,1));
mask_rsum = sum(mask,1);
mask_z = reshape(sum(mask_rsum,2),[k,1]);

% Set indices for trimming
r1 = 70; r2 = 100;
c1 = 63; c2 = 95;
z1 = 1; z2 = 25;

%%% Mask and trim all time steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for q = 1:1:t
    % Get current velocities
    xi = X;
    yi = Y;
    zi = Z;
    ui = Vx(:,:,:,q);
    vi = Vy(:,:,:,q);
    wi = Vz(:,:,:,q);
    
    % Mask the velocities
    um = mask.*ui;
    vm = mask.*vi;
    wm = mask.*wi;
    
    % Zero outlier velocities
    vmag = sqrt(um.^2 + vm.^2 + wm.^2);
    %um(vmag > 0.5) = 0;
    %vm(vmag > 0.5) = 0;
    %wm(vmag > 0.5) = 0;
    
    % Trim the coordinates and velocities
    x = xi(r1:r2,c1:c2,z1:z2);
    y = yi(r1:r2,c1:c2,z1:z2);
    z = zi(r1:r2,c1:c2,z1:z2);
    u = um(r1:r2,c1:c2,z1:z2);
    v = vm(r1:r2,c1:c2,z1:z2);
    w = wm(r1:r2,c1:c2,z1:z2);
    
    % Save the trimmed and masked file as mat file
    save([outpath,out_file,num2str(q,'%05i'),'.mat'],'x','y','z','u','v','w');
    
    % Save output as dat file
    xl = reshape(x,[size(x,1)*size(x,2)*size(x,3),1]);
    yl = reshape(y,[size(x,1)*size(x,2)*size(x,3),1]);
    zl = reshape(z,[size(x,1)*size(x,2)*size(x,3),1]);
    ul = reshape(u,[size(x,1)*size(x,2)*size(x,3),1]);
    vl = reshape(v,[size(x,1)*size(x,2)*size(x,3),1]);
    wl = reshape(w,[size(x,1)*size(x,2)*size(x,3),1]);
    fid = fopen([outpath,out_file,num2str(q,'%05i'),'.dat'],'w');
    fprintf(fid,['TITLE = "',out_file,num2str(q,'%05i'),'"']);
    fprintf(fid,'\nVARIABLES = "X", "Y", "Z", "U", "V", "W"\n');
    %fprintf(fid,['ZONE T="Frame 0", I=',num2str(size(x,1)),', J=',num2str(size(x,2)),', K=',num2str(size(x,3))]);


    % Iterate through all points and write only those that are non-zero
    for zz = 1:1:length(ul)
        cx = xl(zz);
        cy = yl(zz);
        cz = zl(zz);
        cu = ul(zz);
        cv = vl(zz);
        cw = wl(zz);
        fprintf(fid,'\n%.6f %.6f %.6f %.6f %.6f %.6f',cx,cy,cz,cu,cv,cw);
    end
    fclose(fid);
end

