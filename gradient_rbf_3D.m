function [dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz] = gradient_rbf_3D(X,Y,Z,U,V,W,velmask,winsize)
% Inputs should be in physcial coordinates (m and m/s)
% X,Y,Z,U,V,W should be 3D in space matrices 
% velmask - 3D array with NaN's outside mask and 0's inside mask
% winsize - Single value, indicating the size of the window to use. The
% number of points included in the interpolation will be controlled to be
% the cube of the window size

% Initialize the outputs
dudx = zeros(size(X)); dudy = zeros(size(X)); dudz = zeros(size(X));
dvdx = zeros(size(X)); dvdy = zeros(size(X)); dvdz = zeros(size(X));
dwdx = zeros(size(X)); dwdy = zeros(size(X)); dwdz = zeros(size(X));

%%% Iterate through each point in the velocity field
num_iter = sum(sum(sum(~isnan(velmask))));
igrid_init = floor(winsize/2);
num_grid_pts = winsize^3;
iter = 1; prev_mark = 0;
fprintf('Computing velocity gradients...')
for ii = 1:1:size(X,1)
    for jj = 1:1:size(X,2)
        for kk = 1:1:size(X,3)
            x = X(ii,jj,kk);
            y = Y(ii,jj,kk);
            z = Z(ii,jj,kk);
            is_valid = ~isnan(velmask(ii,jj,kk));
            %is_edge = ~isempty(edge_normal{ii,jj,kk});
            if is_valid
                %%% Obtain the neighborhood of points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % The point neighborhood is controlled to always include
                % the same number of points       
                curr_grid_pts = 0;
                igrid = igrid_init;
                while curr_grid_pts < num_grid_pts
                    % Get the indices to take for the current point based
                    % on the window size
                    ind_iiS = max(ii-igrid,1);
                    ind_iiE = min(ii+igrid,size(X,1));
                    ind_jjS = max(jj-igrid,1);
                    ind_jjE = min(jj+igrid,size(X,2));
                    ind_kkS = max(kk-igrid,1);
                    ind_kkE = min(kk+igrid,size(X,3));

                    % Obtain the mesh grid
                    Xmesh = X(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                    Ymesh = Y(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                    Zmesh = Z(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                    Umesh = U(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                    Vmesh = V(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                    Wmesh = W(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                    MASKmesh = velmask(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);

                    % Reshape the mesh grid to a column vector
                    xmesh = reshape(Xmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
                    ymesh = reshape(Ymesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
                    zmesh = reshape(Zmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
                    umesh = reshape(Umesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
                    vmesh = reshape(Vmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
                    wmesh = reshape(Wmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
                    maskmesh = reshape(MASKmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);

                    % Keep only the valid grid points
                    xmesh = xmesh(~isnan(maskmesh));
                    ymesh = ymesh(~isnan(maskmesh));
                    zmesh = zmesh(~isnan(maskmesh));
                    umesh = umesh(~isnan(maskmesh));
                    vmesh = vmesh(~isnan(maskmesh));
                    wmesh = wmesh(~isnan(maskmesh));
                    
                    % Keep only the number of grid points indicated by
                    % the window size to maintain consistency in the
                    % calculation. Keep the closest points.
                    if length(xmesh) > num_grid_pts
                        mesh_dist = sqrt((xmesh-x).^2+(ymesh-y).^2+(zmesh-z).^2); % Distance of each mesh point
                        [~,sinds] = sort(mesh_dist,'ascend'); % Sort the distance from close to far
                        kp_inds = sort(sinds(1:num_grid_pts),'ascend'); % Indices to keep in mesh
                        xmesh = xmesh(kp_inds);
                        ymesh = ymesh(kp_inds);
                        zmesh = zmesh(kp_inds);
                        umesh = umesh(kp_inds);
                        vmesh = vmesh(kp_inds);
                        wmesh = wmesh(kp_inds);
                    end

                    % Increment the grid step out value
                    igrid = igrid + 1;
                    % Recompute the number of grid points
                    curr_grid_pts = length(xmesh);
                end
% % %                 else
% % %                     % If the point is an edge, determine the points in the
% % %                     % grid based on the unit normal
% % %                     unorm = abs(edge_normal{ii,jj,kk});
% % %                     % Set igrid, jgrid, and kgrid values based on the
% % %                     % normal direction (want to extend most into the middle
% % %                     % of the flow field)
% % %                     i_pts = round(unorm(1) * num_grid_pts);
% % %                     j_pts = round(unorm(2) * num_grid_pts);
% % %                     k_pts = round(unorm(3) * num_grid_pts);
% % %                     % Set the i, j, k window sizes
% % %                     iwin = 2;%max(i_pts/5,2);
% % %                     jwin = 2;%max(j_pts/5,2);
% % %                     kwin = 2;%max(k_pts/5,2);
% % %                     
% % %                     % Get the indices to take for the current point based
% % %                     % on the window size
% % %                     ind_iiS = max(ii-floor(iwin/2),1);
% % %                     ind_iiE = min(ii+ceil(iwin/2),size(X,1));
% % %                     ind_jjS = max(jj-floor(jwin/2),1);
% % %                     ind_jjE = min(jj+ceil(jwin/2),size(X,2));
% % %                     ind_kkS = max(kk-floor(kwin/2),1);
% % %                     ind_kkE = min(kk+ceil(kwin/2),size(X,3));
% % % 
% % %                     % Obtain the mesh grid
% % %                     Xmesh = X(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % %                     Ymesh = Y(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % %                     Zmesh = Z(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % %                     Umesh = U(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % %                     Vmesh = V(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % %                     Wmesh = W(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % %                     MASKmesh = velmask(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
% % % 
% % %                     % Reshape the mesh grid to a column vector
% % %                     xmesh = reshape(Xmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % %                     ymesh = reshape(Ymesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % %                     zmesh = reshape(Zmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % %                     umesh = reshape(Umesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % %                     vmesh = reshape(Vmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % %                     wmesh = reshape(Wmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % %                     maskmesh = reshape(MASKmesh,[size(Xmesh,1)*size(Xmesh,2)*size(Xmesh,3),1]);
% % % 
% % %                     % Keep only the valid grid points
% % %                     xmesh = xmesh(~isnan(maskmesh));
% % %                     ymesh = ymesh(~isnan(maskmesh));
% % %                     zmesh = zmesh(~isnan(maskmesh));
% % %                     umesh = umesh(~isnan(maskmesh));
% % %                     vmesh = vmesh(~isnan(maskmesh));
% % %                     wmesh = wmesh(~isnan(maskmesh));
% % %                 end
% % %             
                %%% Compute the gradients using thin plate spline (RBF-TPS) %%%%%%%%%%%
                [c_dudx,c_dudy,c_dudz,c_dvdx,c_dvdy,c_dvdz,c_dwdx,c_dwdy,c_dwdz] = gradient_tps3_3D(x,y,z,xmesh,ymesh,zmesh,umesh,vmesh,wmesh);
                
                %%% Save the gradient in the matrix
                dudx(ii,jj,kk) = c_dudx;
                dudy(ii,jj,kk) = c_dudy;
                dudz(ii,jj,kk) = c_dudz;
                dvdx(ii,jj,kk) = c_dvdx;
                dvdy(ii,jj,kk) = c_dvdy;
                dvdz(ii,jj,kk) = c_dvdz;
                dwdx(ii,jj,kk) = c_dwdx;
                dwdy(ii,jj,kk) = c_dwdy;
                dwdz(ii,jj,kk) = c_dwdz;
                
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
fprintf(' complete \n')

end % End of function

                    
