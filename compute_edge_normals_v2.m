function [edge_normal] = compute_edge_normals_v2(velmask,xmask,ymask,zmask,stlmask,Xstl,Ystl,Zstl,SMOOTHNORMS,PLOTNORMS)
% Velocity mask - 3D matrix containing 0's for all points within the mask
% and NaN's for all points outside of the mask
% x_mask - 1D array containing the dimensions of the x-coord (cols) of the mask
% y_mask - 1D array containing the dimensions of the y-coord (rows) of the mask
% z_mask - 1D array containing the dimensions of the z-coord (planes) of the mask
% stlmask - 3D matrix containing the STL mask, of the same format as the
% velocity mask
% SMOOTHNORMS - Choose whether or not to smooth the normals (using 2-pass
% UOD), always suggested to be on
% PLOTNORMS - Choose whether or not to plot the normals as they are computed, 
% (this is mostly for debugging. If it is off, the final normals will be
% plotted to the user still.)

%%% Find the edges of the velocity grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xm,Ym,Zm] = meshgrid(xmask,ymask,zmask);
edge_pts = zeros(size(Zm));
for ii = 2:1:size(Xm,1)-1
    for jj = 2:1:size(Xm,2)-1
        for kk = 2:1:size(Xm,3)-1
            % Determine if the current point is in the boundary
            if ~isnan(velmask(ii,jj,kk))
                % Determine if point is an edge
                x_neighs = ~isnan(velmask(ii-1,jj,kk)) + ~isnan(velmask(ii+1,jj,kk));     
                y_neighs = ~isnan(velmask(ii,jj-1,kk)) + ~isnan(velmask(ii,jj+1,kk));
                z_neighs = ~isnan(velmask(ii-1,jj,kk-1)) + ~isnan(velmask(ii,jj,kk+1));
                all_neighs_gd = x_neighs + y_neighs + z_neighs;
                if all_neighs_gd < 6
                    % Set the edge point to 1
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


%%% Find the edges of the STL grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xstl = permute(Xstl(1,:,1),[2 1 3]);
ystl = Ystl(:,1,1);
zstl = permute(Zstl(1,1,:),[3,1,2]);
stl_edge_pts = zeros(size(Zstl));
for ii = 2:1:size(Xstl,1)-1
    for jj = 2:1:size(Xstl,2)-1
        for kk = 2:1:size(Xstl,3)-1
            % Determine if the current point is in the boundary
            if ~isnan(stlmask(ii,jj,kk))
                % Determine if point is an edge
                x_neighs = ~isnan(stlmask(ii-1,jj,kk)) + ~isnan(stlmask(ii+1,jj,kk));     
                y_neighs = ~isnan(stlmask(ii,jj-1,kk)) + ~isnan(stlmask(ii,jj+1,kk));
                z_neighs = ~isnan(stlmask(ii-1,jj,kk-1)) + ~isnan(stlmask(ii,jj,kk+1));
                all_neighs_gd = x_neighs + y_neighs + z_neighs;
                if all_neighs_gd < 6
                    % Set the edge point to 1
                    stl_edge_pts(ii,jj,kk) = 1;
                end                        
            end
        end
    end
end
% Obtain all boundary points as a line -  for plotting purposes
stl_x_all_pts = reshape(Xstl,[size(Xstl,1)*size(Xstl,2)*size(Xstl,3),1]);
stl_y_all_pts = reshape(Ystl,[size(Ystl,1)*size(Ystl,2)*size(Ystl,3),1]);
stl_z_all_pts = reshape(Zstl,[size(Zstl,1)*size(Zstl,2)*size(Zstl,3),1]);
stl_is_edge = reshape(stl_edge_pts,[size(Xstl,1)*size(Xstl,2)*size(Xstl,3),1]);
% Limit to the points in the mask only
stl_x_epts = stl_x_all_pts(stl_is_edge == 1);%(~isnan(mask_pts));
stl_y_epts = stl_y_all_pts(stl_is_edge == 1); %~isnan(mask_pts));
stl_z_epts = stl_z_all_pts(stl_is_edge == 1); %~isnan(mask_pts));

%%% Compute the wall normals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each boundary point, compute its normal vector
% This is done using the STL by finding the closest point on the STL and
% the associated STL neighborhood
fprintf('\nComputing wall normals...')
edge_normal = cell(size(edge_pts));
stencil_size = 5;
half_stencil = floor(stencil_size/2);
num_neigh_pts = 25;
% Plot for debugging
if PLOTNORMS, figure(27); hold off; scatter3(x_epts,y_epts,z_epts); hold on; end
% Line t array
line_t = linspace(0.5,10,50);
neigh_t = linspace(0,5,20);
iter = 1; prev_mark = 0;
num_iter = length(x_epts);
for ii = 1:1:size(Xm,1)
    for jj = 1:1:size(Xm,2)
        for kk = 1:1:size(Xm,3)
            % Only calculate if the point is an edge
            if edge_pts(ii,jj,kk) == 1
                %%% Get current edge point info %%%%%%%%%%%%%%%%%%%
                x = Xm(ii,jj,kk);
                y = Ym(ii,jj,kk);
                z = Zm(ii,jj,kk);
                
                if ii == 14 && jj == 12 && kk == 3
                    keyboard
                end
                
                %%% Get the 50 closest boundary points %%%%%%%%%%%%%%%%%%%%
                % The STL points are used for this step - larger number of
                % points are used because the STL has high resolution
                % Compute the point distance
                pt_dists = sqrt((stl_x_epts-x).^2+(stl_y_epts-y).^2+(stl_z_epts-z).^2);
                [srt_dists,sinds] = sort(pt_dists,'ascend');
                % Add the points through an iterative process to
                % pick the points closest to the panel of points
                kpX = stl_x_epts(sinds(1:2));
                kpY = stl_y_epts(sinds(1:2));
                kpZ = stl_z_epts(sinds(1:2));
                kept_inds = sinds(1:2);
                pt_kept = zeros(size(stl_x_epts));
                pt_kept(sinds(1:2)) = 1;
                max_conn_grid_dist = sqrt((xstl(2)-xstl(1))^2 + (ystl(2)-ystl(1))^2 + (zstl(2)-zstl(1))^2);
                for qq = 3:1:num_neigh_pts
                    all_pt_dists = 100*ones(size(stl_x_epts));
                    for zz = 1:1:length(kpX)
                        all_dists = sqrt((stl_x_epts-kpX(zz)).^2+(stl_y_epts-kpY(zz)).^2+(stl_z_epts-kpZ(zz)).^2);
                        all_dists(pt_kept == 1) = 100;
                        all_dists(all_dists > max_conn_grid_dist) = 5*all_dists(all_dists > max_conn_grid_dist);
                        all_dists = all_dists + pt_dists;
                        all_pt_dists = min([all_dists,all_pt_dists],[],2);
                    end
                    [srt_dists,sinds] = sort(all_pt_dists,'ascend');
                    % Add new point
                    kpX = [kpX;stl_x_epts(sinds(1))];
                    kpY = [kpY;stl_y_epts(sinds(1))];
                    kpZ = [kpZ;stl_z_epts(sinds(1))];
                    kept_inds = [kept_inds;sinds(1)];
                    pt_kept(sinds(1)) = 1;
                end
                ptsurfX = kpX;
                ptsurfY = kpY;
                ptsurfZ = kpZ;

                % Plot for debugging
                if PLOTNORMS
                    figure(26); hold off; scatter3(stl_x_epts,stl_y_epts,stl_z_epts,2,'MarkerEdgeColor',[0.6,0.6,0.6])%,'MarkerFaceColor',[0.6,0.6,0.6])
                    hold on; scatter3(x_epts,y_epts,z_epts,8,'b','MarkerFaceColor','b')
                    scatter3(kpX,kpY,kpZ,8,'r','MarkerFaceColor','r'); 
                    scatter3(x,y,z,10,'k','MarkerFaceColor','k'); drawnow;
                end

                %%% Compute the normal values %%%%%%%%%%%%%%%%%%%%%
                % Determine if the mesh should be X/Y, Y/Z, or X/Z
                % based. This is chosen based on how many points
                % exist in each plane. The maximum number of points
                % is selected.
                % Number of XY points
                mat1 = [ptsurfX,ptsurfY];
                unqmat = unique(mat1,'rows');
                numXY = size(unqmat,1);
                % Number of YZ points
                mat1 = [ptsurfY,ptsurfZ];
                unqmat = unique(mat1,'rows');
                numYZ = size(unqmat,1);
                mat1 = [ptsurfX,ptsurfZ];
                unqmat = unique(mat1,'rows');
                numXZ = size(unqmat,1);
                [~,opt_plane] = max([numXY,numYZ,numXZ]);
                xmesh = linspace(min(ptsurfX),max(ptsurfX),25); %unique(ptsurfX); %sort(ptsurfX,'ascend'); %
                ymesh = linspace(min(ptsurfY),max(ptsurfY),25); %unique(ptsurfY); %sort(ptsurfY,'ascend');
                zmesh = linspace(min(ptsurfZ),max(ptsurfZ),25); %unique(ptsurfZ); %sort(ptsurfZ,'ascend'); %
                ptsurfX_adj = 0; ptsurfY_adj = 0; ptsurfZ_adj = 0;
                if opt_plane == 1
                    % Use an X/Y based plane
                    [Xmesh,Ymesh] = meshgrid(xmesh,ymesh);
                    % Adjust the surface points so there are no
                    % duplicates
                    all_pts_inc = zeros(size(ptsurfX));
                    ct = 1; cind = 1;
                    while sum(all_pts_inc==0) > 0
                        % Get current point
                        currX = ptsurfX(cind);
                        currY = ptsurfY(cind);
                        % Get instances of current point
                        curr_indX = find(ptsurfX == currX);
                        curr_indY = find(ptsurfY == currY);
                        curr_inds = intersect(curr_indX,curr_indY);
                        % Add the current index/indices to the
                        % adjusted points list
                        ptsurfX_adj(ct,1) = currX;
                        ptsurfY_adj(ct,1) = currY;
                        ptsurfZ_adj(ct,1) = mean(ptsurfZ(curr_inds));

                        % Add the current point to the list of
                        % included points
                        all_pts_inc(curr_inds) = 1;

                        % Update the current point and counter
                        cind = find(all_pts_inc == 0,1,'first');
                        ct = ct + 1;
                    end
                    % Interpolate to the grid
                    Zmesh = griddata(ptsurfX_adj,ptsurfY_adj,ptsurfZ_adj,Xmesh,Ymesh,'cubic');
                    % Get the outward normals for the surface
                    [nx,ny,nz] = surfnorm(Xmesh,Ymesh,Zmesh);
                    % Reshape the normal and mesh arrays
                    Nx = reshape(nx,[size(nx,1)*size(nx,2),1]);
                    Ny = reshape(ny,[size(ny,1)*size(ny,2),1]);
                    Nz = reshape(nz,[size(nz,1)*size(nz,2),1]);
                    xml = reshape(Xmesh,[size(Xmesh,1)*size(Xmesh,2),1]);
                    yml = reshape(Ymesh,[size(Ymesh,1)*size(Ymesh,2),1]);
                    zml = reshape(Zmesh,[size(Zmesh,1)*size(Zmesh,2),1]);
                    % Keep only the non NaN's
                    kp_inds = ~isnan(Nx);
                    Nx = Nx(kp_inds);
                    Ny = Ny(kp_inds);
                    Nz = Nz(kp_inds);
                    xml = xml(kp_inds);
                    yml = yml(kp_inds);
                    zml = zml(kp_inds);
                    % Find the normal value closest to the point of interest
                    [norm_dist,opt_ind] = min(sqrt((xml-x).^2 + (yml-y).^2));
                    % Get the closest normal, negate to make inward
                    normX = -Nx(opt_ind);
                    normY = -Ny(opt_ind);
                    normZ = -Nz(opt_ind);
                elseif opt_plane == 2
                    % Use a Y/Z based plane
                    [Ymesh,Zmesh] = meshgrid(ymesh,zmesh);
                    % Adjust the surface points so there are no
                    % duplicates
                    all_pts_inc = zeros(size(ptsurfX));
                    ct = 1; cind = 1;
                    while sum(all_pts_inc==0) > 0
                        % Get current point
                        currY = ptsurfY(cind);
                        currZ = ptsurfZ(cind);
                        % Get instances of current point
                        curr_indY = find(ptsurfY == currY);
                        curr_indZ = find(ptsurfZ == currZ);
                        curr_inds = intersect(curr_indY,curr_indZ);
                        % Add the current index/indices to the
                        % adjusted points list
                        ptsurfY_adj(ct,1) = currY;
                        ptsurfZ_adj(ct,1) = currZ;
                        ptsurfX_adj(ct,1) = mean(ptsurfX(curr_inds));

                        % Add the current point to the list of
                        % included points
                        all_pts_inc(curr_inds) = 1;

                        % Update the current point and counter
                        cind = find(all_pts_inc == 0,1,'first');
                        ct = ct + 1;
                    end

                    % Interpolate to the grid
                    Xmesh = griddata(ptsurfY_adj,ptsurfZ_adj,ptsurfX_adj,Ymesh,Zmesh,'cubic');
                    % Get the outward normals for the surface
                    [nx,ny,nz] = surfnorm(Xmesh,Ymesh,Zmesh);
                    % Reshape the normal and mesh arrays
                    Nx = reshape(nx,[size(nx,1)*size(nx,2),1]);
                    Ny = reshape(ny,[size(ny,1)*size(ny,2),1]);
                    Nz = reshape(nz,[size(nz,1)*size(nz,2),1]);
                    xml = reshape(Xmesh,[size(Xmesh,1)*size(Xmesh,2),1]);
                    yml = reshape(Ymesh,[size(Ymesh,1)*size(Ymesh,2),1]);
                    zml = reshape(Zmesh,[size(Zmesh,1)*size(Zmesh,2),1]);
                    % Keep only the non NaN's
                    kp_inds = ~isnan(Nx);
                    Nx = Nx(kp_inds);
                    Ny = Ny(kp_inds);
                    Nz = Nz(kp_inds);
                    xml = xml(kp_inds);
                    yml = yml(kp_inds);
                    zml = zml(kp_inds);
                    % Find the normal value closest to the point of interest
                    [norm_dist,opt_ind] = min(sqrt((yml-y).^2 + (zml-z).^2));
                    % Get the closest normal, negate to make inward
                    normX = -Nx(opt_ind);
                    normY = -Ny(opt_ind);
                    normZ = -Nz(opt_ind);
                else
                    % Use an X/Z based plane
                    [Xmesh,Zmesh] = meshgrid(xmesh,zmesh);
                    % Adjust the surface points so there are no
                    % duplicates
                    all_pts_inc = zeros(size(ptsurfX));
                    ct = 1; cind = 1;
                    while sum(all_pts_inc==0) > 0
                        % Get current point
                        currX = ptsurfX(cind);
                        currZ = ptsurfZ(cind);
                        % Get instances of current point
                        curr_indX = find(ptsurfX == currX);
                        curr_indZ = find(ptsurfZ == currZ);
                        curr_inds = intersect(curr_indX,curr_indZ);
                        % Add the current index/indices to the
                        % adjusted points list
                        ptsurfX_adj(ct,1) = currX;
                        ptsurfZ_adj(ct,1) = currZ;
                        ptsurfY_adj(ct,1) = mean(ptsurfY(curr_inds));

                        % Add the current point to the list of
                        % included points
                        all_pts_inc(curr_inds) = 1;

                        % Update the current point and counter
                        cind = find(all_pts_inc == 0,1,'first');
                        ct = ct + 1;
                    end

                    % Interpolate to the grid
                    Ymesh = griddata(ptsurfX_adj,ptsurfZ_adj,ptsurfY_adj,Xmesh,Zmesh,'cubic');
                    % Get the outward normals for the surface
                    [nx,ny,nz] = surfnorm(Xmesh,Ymesh,Zmesh);
                    % Reshape the normal and mesh arrays
                    Nx = reshape(nx,[size(nx,1)*size(nx,2),1]);
                    Ny = reshape(ny,[size(ny,1)*size(ny,2),1]);
                    Nz = reshape(nz,[size(nz,1)*size(nz,2),1]);
                    xml = reshape(Xmesh,[size(Xmesh,1)*size(Xmesh,2),1]);
                    yml = reshape(Ymesh,[size(Ymesh,1)*size(Ymesh,2),1]);
                    zml = reshape(Zmesh,[size(Zmesh,1)*size(Zmesh,2),1]);
                    % Keep only the non NaN's
                    kp_inds = ~isnan(Nx);
                    Nx = Nx(kp_inds);
                    Ny = Ny(kp_inds);
                    Nz = Nz(kp_inds);
                    xml = xml(kp_inds);
                    yml = yml(kp_inds);
                    zml = zml(kp_inds);
                    % Find the normal value closest to the point of interest
                    [norm_dist,opt_ind] = min(sqrt((xml-x).^2 + (zml-z).^2));
                    % Get the closest normal, negate to make inward
                    normX = -Nx(opt_ind);
                    normY = -Ny(opt_ind);
                    normZ = -Nz(opt_ind);
                end

                %%% Get normal direction %%%%%%%%%%%%%%%%%%%%%%%%%%
                % Determine if the direction of the unit normal is
                % correct, or if it needs to be flipped. The normal
                % should not point into the masked NaN values
                % Get line of points going in direction of the
                % normal, and the direction opposite the normal
                line_forw_x = x + normX*line_t;
                line_forw_y = y + normY*line_t;
                line_forw_z = z + normZ*line_t;
                line_back_x = x - normX*line_t;
                line_back_y = y - normY*line_t;
                line_back_z = z - normZ*line_t;
                % Convert the line points to the nearest pixel values
                forw_inds = [];
                back_inds = [];
                for ptnum = 1:1:length(line_forw_x)
                    % Current line point
                    clx = line_forw_x(ptnum);
                    cly = line_forw_y(ptnum);
                    clz = line_forw_z(ptnum);
                    % Nearest STL point
                    if (clx>=min(xstl)) && (clx<=max(xstl)) && (cly>=min(ystl)) && (cly<=max(ystl)) && (clz>=min(zstl)) && (clz<=max(zstl))
                        [~,mx_ind] = min(abs(xstl-clx));
                        [~,my_ind] = min(abs(ystl-cly));
                        [~,mz_ind] = min(abs(zstl-clz));
                        forw_inds = [forw_inds;my_ind,mx_ind,mz_ind]; % Order as (row,col,z)
                    end
                    % Current line point
                    clx = line_back_x(ptnum);
                    cly = line_back_y(ptnum);
                    clz = line_back_z(ptnum);
                    % Nearest STL point
                    if (clx>=min(xstl)) && (clx<=max(xstl)) && (cly>=min(ystl)) && (cly<=max(ystl)) && (clz>=min(zstl)) && (clz<=max(zstl))
                        [~,mx_ind] = min(abs(xstl-clx));
                        [~,my_ind] = min(abs(ystl-cly));
                        [~,mz_ind] = min(abs(zstl-clz));
                        back_inds = [back_inds;my_ind,mx_ind,mz_ind]; % Order as (row,col,z)
                    end
                end
                % Keep only the unique rows of points
                [~,forw_rows] = unique(forw_inds,'rows');
                [~,back_rows] = unique(back_inds,'rows');
                forw_rows = sort(forw_rows,'ascend');
                back_rows = sort(back_rows,'ascend');
                forw_inds = forw_inds(forw_rows,:);
                back_inds = back_inds(back_rows,:);
                forw_pixline = zeros(size(forw_inds,1),1);
                back_pixline = zeros(size(back_inds,1),1);
                % Get the forward and backward pixel values
                for rownum = 1:1:size(forw_inds,1)
                    forw_pixline(rownum) = stlmask(forw_inds(rownum,1),forw_inds(rownum,2),forw_inds(rownum,3));
                end
                for rownum = 1:1:size(back_inds,1)
                    back_pixline(rownum) = stlmask(back_inds(rownum,1),back_inds(rownum,2),back_inds(rownum,3));
                end
                forw_nan = find(isnan(forw_pixline),1,'first');
                if isempty(forw_nan),forw_nan = length(back_pixline)+1; end
                back_nan = find(isnan(back_pixline),1,'first');
                if isempty(back_nan),back_nan = length(forw_pixline)+1; end

                % Determine if the normal needs to be flipped
                if back_nan > forw_nan
                    normX = -normX;
                    normY = -normY;
                    normZ = -normZ;
                end

                %%% Save the unit normal %%%%%%%%%%%%%%%%%%%%%%%%%%
                unorm = [normX,normY,normZ];
                unorm = unorm/sum(sqrt(unorm.^2));
                edge_normal{ii,jj,kk} = unorm;

                % Plot for debugging
                if PLOTNORMS
                    figure(27); hold on; quiver3(x,y,z,unorm(1),unorm(2),unorm(3),'r'); drawnow
                end

                % Print the percent completed for steps of 10
                perc_complete = floor(100*iter/num_iter);
                if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                    fprintf(' %i%%',perc_complete)
                    prev_mark = perc_complete;
                end

                % Increment the iteration number
                iter = iter + 1;
            end
        end
    end
end

if SMOOTHNORMS
    fprintf('\nRunning outlier detection smoothing on normals...')
    % Run a UOD on all normals to determine if there are any erroneous
    % normals computed
    max_conn_grid_dist = sqrt((xmask(2)-xmask(1))^2 + (ymask(2)-ymask(1))^2 + (zmask(2)-zmask(1))^2);
    grid_res = abs(xmask(2)-xmask(1));
    num_passes = 2;
    pass_tol = [3,2];
    e = 0.1*grid_res;
    iter = 1; prev_mark = 0;
    num_iter = num_passes*length(x_epts);
    eval = zeros(size(x_epts));
    if PLOTNORMS, figure(29); hold off; scatter3(x_epts,y_epts,z_epts); end
    for pass_num = 1:1:num_passes
        tol = pass_tol(pass_num);
        for ii = 1:1:size(Xm,1)
            for jj = 1:1:size(Xm,2)
                for kk = 1:1:size(Xm,3)
                    %%% Get current edge point info %%%%%%%%%%%%%%%%%%%
                    x = Xm(ii,jj,kk);
                    y = Ym(ii,jj,kk);
                    z = Zm(ii,jj,kk);
                    curr_normal = edge_normal{ii,jj,kk};
                    
                    if ~isempty(curr_normal)
                        %%% Find the 15 closest normal points %%%%%%%%%%%%%
                        % Compute the point distance
                        pt_dists = sqrt((x_epts-x).^2+(y_epts-y).^2+(z_epts-z).^2);
                        [srt_dists,sinds] = sort(pt_dists,'ascend');
                        % Add the points through an iterative process to
                        % pick the points closest to the panel of points
                        kpX = x_epts(sinds(1:2));
                        kpY = y_epts(sinds(1:2));
                        kpZ = z_epts(sinds(1:2));
                        kept_inds = sinds(1:2);
                        pt_kept = zeros(size(x_epts));
                        pt_kept(sinds(1:2)) = 1;
                        num_neighs = 15;
                        for qq = 3:1:num_neighs
                            all_pt_dists = 100*ones(size(x_epts));
                            for zz = 1:1:length(kpX)
                                all_dists = sqrt((x_epts-kpX(zz)).^2+(y_epts-kpY(zz)).^2+(z_epts-kpZ(zz)).^2);
                                all_dists(pt_kept == 1) = 100;
                                all_dists(all_dists > max_conn_grid_dist) = 5*all_dists(all_dists > max_conn_grid_dist);
                                all_dists = all_dists + pt_dists;
                                all_pt_dists = min([all_dists,all_pt_dists],[],2);
                            end
                            [srt_dists,sinds] = sort(all_pt_dists,'ascend');
                            % Add new point
                            kpX = [kpX;x_epts(sinds(1))];
                            kpY = [kpY;y_epts(sinds(1))];
                            kpZ = [kpZ;z_epts(sinds(1))];
                            kept_inds = [kept_inds;sinds(1)];

                            pt_kept(sinds(1)) = 1;
                        end

                        %%% Get the normal for each neighbor point %%%%%%%%
                        neigh_norms = zeros(num_neighs,3);
                        neigh_norms(1,:) = curr_normal;
                        if PLOTNORMS
                           figure(28); hold off; scatter3(kpX,kpY,kpZ); hold on; scatter3(x,y,z);
                           xlabel('X'); ylabel('Y'); zlabel('Z');
                           quiver3(x,y,z,neigh_norms(1,1),neigh_norms(1,2),neigh_norms(1,3),'r'); drawnow
                        end
                        for qq = 2:1:num_neighs
                            % Current neighbor point
                            cx = kpX(qq);
                            cy = kpY(qq);
                            cz = kpZ(qq);
                            % Find the index in the X,Y,Z
                            [~,x_ind] = min(abs(xmask-cx));
                            [~,y_ind] = min(abs(ymask-cy));
                            [~,z_ind] = min(abs(zmask-cz));
                            % Get the edge normal
                            new_norm = edge_normal{y_ind,x_ind,z_ind};
                            % Add it to the neighbor list
                            neigh_norms(qq,:) = new_norm;
                            if PLOTNORMS
                               figure(28); quiver3(cx,cy,cz,new_norm(1),new_norm(2),new_norm(3),'k'); drawnow
                            end
                        end        

                        %%% Compute the median normal %%%%%%%%%%%%%%%%%%%%%
                        med_norm = median(neigh_norms,1);
                        med_norm_mag = sqrt(sum(med_norm.^2));
                        med_norm = med_norm/med_norm_mag; % Make the median normal a unit normal
                        % Compute the difference between the normal and the
                        % median normal
                        norm_error = sum(abs(med_norm - curr_normal));
                        % Compute residual and redisual median
                        ures = abs(neigh_norms(:,1) - med_norm(1));
                        vres = abs(neigh_norms(:,2) - med_norm(2));
                        wres = abs(neigh_norms(:,3) - med_norm(3));
                        rumed = nanmedian(ures(:));
                        rvmed = nanmedian(vres(:));
                        rwmed = nanmedian(wres(:));

                        % Compute residual of current point
                        ru = abs(curr_normal(1) - med_norm(1))/(rumed + e);
                        rv = abs(curr_normal(2) - med_norm(2))/(rvmed + e);
                        rw = abs(curr_normal(3) - med_norm(3))/(rwmed + e);
                        restot = sqrt(ru.^2 + rv.^2 + rw.^2); % Total residual

                        %%% Run UOD style smoothing %%%%%%%%%%%%%%%%%%%%%%%
                        % Change to median normal if the total residual is
                        % greater than the tolerance
                        if restot > tol
                            edge_normal{ii,jj,kk} = med_norm;
                            eval(kept_inds(1)) = pass_num;
                            curr_normal = med_norm;
                        end

                        % Plot for debugging
                        if PLOTNORMS
                            if pass_num == 1
                                figure(29); hold on; quiver3(x,y,z,curr_normal(1),curr_normal(2),curr_normal(3),'r'); drawnow;
                            else
                                figure(29); hold on; quiver3(x,y,z,curr_normal(1),curr_normal(2),curr_normal(3),'k'); drawnow;
                            end
                            figure(28); hold on; quiver3(x,y,z,curr_normal(1),curr_normal(2),curr_normal(3),'g');
                            title(['residual = ',num2str(restot,'%.2f')]); drawnow;
                        end

                        % Print the percent completed for steps of 10
                        perc_complete = floor(100*iter/num_iter);
                        if (mod(perc_complete,10) == 0) && (perc_complete ~= prev_mark)
                            fprintf(' %i%%',perc_complete)
                            prev_mark = perc_complete;
                        end

                        % Increment the iteration
                        iter = iter + 1;
                    end
                end
            end
        end
    end
end

fprintf('\n')
end

