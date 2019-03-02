function [wss_mag,wss_x,wss_y,wss_z] = compute_wss(edge_normal,viscosity_mu,grid_size,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
% Computes the 3D wall shear stress from the input normals and gradients
% This computation is described in the Appendix of Stalder et al. (2008)
% grid_size is used to take a median (or max) derivative value in case the
% boundary lowers this derivative - if it is observed that the gradient is
% losing magnitude near the wall, this should be set to reduce that effect.
% A grid_size of 0 corresponds to looking at a point

% Initialize the outputs
wss_mag = zeros(size(dudx));
wss_x = zeros(size(dudx));
wss_y = zeros(size(dudx));
wss_z = zeros(size(dudx));

% Set the grid window based on the grid_size input
gwin = floor(grid_size/2);

% Iterate through all points
for ii = 1:1:size(dudx,1)
    for jj = 1:1:size(dudx,2)
        for kk = 1:1:size(dudx,3)
            is_edge = ~isempty(edge_normal{ii,jj,kk});
            
            %if ii == 83 && jj == 60
            %    keyboard
            %end
            
            % Only continue if the point is an edge
            if is_edge
                % Get the gradients at the current point - This is done
                % through a window process
                % Get the window indices
                ind_iiS = max(ii-gwin,1);
                ind_iiE = min(ii+gwin,size(dudx,1));
                ind_jjS = max(jj-gwin,1);
                ind_jjE = min(jj+gwin,size(dudx,2));
                ind_kkS = max(kk-gwin,1);
                ind_kkE = min(kk+gwin,size(dudx,3));
                % Get the derivative windows
                pt_dudx = dudx(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE); 
                pt_dudy = dudy(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE); 
                pt_dudz = dudz(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                pt_dvdx = dvdx(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE); 
                pt_dvdy = dvdy(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE); 
                pt_dvdz = dvdz(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                pt_dwdx = dwdx(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                pt_dwdy = dwdy(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                pt_dwdz = dwdz(ind_iiS:1:ind_iiE,ind_jjS:1:ind_jjE,ind_kkS:1:ind_kkE);
                % Reshape the window to a line
                pt_dudx = pt_dudx(:);
                pt_dudy = pt_dudy(:);
                pt_dudz = pt_dudz(:);
                pt_dvdx = pt_dvdx(:);
                pt_dvdy = pt_dvdy(:);
                pt_dvdz = pt_dvdz(:);
                pt_dwdx = pt_dwdx(:);
                pt_dwdy = pt_dwdy(:);
                pt_dwdz = pt_dwdz(:);
                % Take the maximum value
                [~,ikp] = max(abs(pt_dudx)); pt_dudx = pt_dudx(ikp);
                [~,ikp] = max(abs(pt_dudy)); pt_dudy = pt_dudy(ikp);
                [~,ikp] = max(abs(pt_dudz)); pt_dudz = pt_dudz(ikp);
                [~,ikp] = max(abs(pt_dvdx)); pt_dvdx = pt_dvdx(ikp);
                [~,ikp] = max(abs(pt_dvdy)); pt_dvdy = pt_dvdy(ikp);
                [~,ikp] = max(abs(pt_dvdz)); pt_dvdz = pt_dvdz(ikp);
                [~,ikp] = max(abs(pt_dwdx)); pt_dwdx = pt_dwdx(ikp);
                [~,ikp] = max(abs(pt_dwdy)); pt_dwdy = pt_dwdy(ikp);
                [~,ikp] = max(abs(pt_dwdz)); pt_dwdz = pt_dwdz(ikp);
                
                % Get the current normal value
                norm = edge_normal{ii,jj,kk};
                
                % Compute the wall shear stress vector
                tau_x = 2*norm(1)*pt_dudx + norm(2)*(pt_dudy+pt_dvdx) + norm(3)*(pt_dudz+pt_dwdx);
                tau_y = norm(1)*(pt_dudy+pt_dvdx) + 2*norm(2)*(pt_dvdy) + norm(3)*(pt_dvdz+pt_dwdy);
                tau_z = norm(1)*(pt_dudz+pt_dwdx) + norm(2)*(pt_dvdz+pt_dwdy) + 2*norm(3)*(pt_dwdz);
                % Compute the wall shear stress magnitude
                tau_mag = sqrt(tau_x^2 + tau_y^2 + tau_z^2);
                
                % Save the wall shear stress into the output matrices
                wss_mag(ii,jj,kk) = viscosity_mu*tau_mag;
                wss_x(ii,jj,kk) = viscosity_mu*tau_x;
                wss_y(ii,jj,kk) = viscosity_mu*tau_y;
                wss_z(ii,jj,kk) = viscosity_mu*tau_z;
            end
        end
    end
end

end

