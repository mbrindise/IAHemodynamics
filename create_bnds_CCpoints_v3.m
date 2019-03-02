%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a mask image using the inputed points and segment.
% Written by: Melissa Brindise
% Date created: November 11, 2017
% --- Inputs --------------------------------------------------------------
% x - the x coordinates of the STL points (should be in single z-plane)
% y - the y coordinates of the STL points (should be in single z-plane)
% pt_segs - the array indicated which connected segment the point belongs
%           to. Must be the same size as x and y.
% cseg - the current segment to evaluate the boundary/mask for
% xp - the x coordinates of the resulting mask created
% yp - the y coordinates of the resulting mask created
% iter - the iteration number of running the mask identification. If iter
%        is greater than 1, the code will have the user select the center point.
% --- OUTPUTS -------------------------------------------------------------
% masked_plane - the mask image created using the inputed points and segment
% xcnt - the x center location used for the polar unwrapping
% ycnt - the y center location used for the polar unwrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [masked_plane] = create_bnds_CCpoints_v3(x,y,pt_segs,cseg,xp,yp,iter,consol_thresh,neigh_thresh)

% Get mask information
[Xm,Ym] = meshgrid(xp,yp);
xml = reshape(Xm,[size(Xm,1)*size(Xm,2),1]);
yml = reshape(Ym,[size(Ym,1)*size(Ym,2),1]);

% Get unique x and y values to create indices
xr = unique(x); % Real coordinate x
yr = unique(y); % Real coordinate y

% Determine which segment has the most number of points in it
num_segs = length(unique(pt_segs)); % Number of segments within the connected points
for q = 1:1:num_segs
    num_pts(q) = sum(pt_segs == q);
end
[~,q] = max(num_pts); % Set the segment to evaluate as the one with the max number of points

% Initialize counters and thresholds
plane_num = 1;
thresh = 1;

% Order the points in the segment
segi = find(pt_segs == cseg);    

% Get the points within the segment only
xs = x(segi);
ys = y(segi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Consolidate points that are too close together %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_thresh = consol_thresh;
pt_consol = zeros(size(xs));
ALL_PTS_CONSOL = 0;
curr_ind = 1;
add_ct = 1;
while ALL_PTS_CONSOL == 0
    % Get current x and y
    curr_x = xs(curr_ind);
    curr_y = ys(curr_ind);
    
    % Compute distances
    dist_inds = sqrt((xs-curr_x).^2 + (ys-curr_y).^2);
    dist_inds = dist_inds + 500*pt_consol; % Add large value for points already consolidated
    
    % Get list of indices within distance threshold
    in_thresh = find(dist_inds <= dist_thresh);
    
    % Consolidate points that are too close to each other
    if length(in_thresh) > 1 % Too many points too close together
        xs_c(add_ct,1) = mean(xs(in_thresh));
        ys_c(add_ct,1) = mean(ys(in_thresh));
    else
        xs_c(add_ct,1) = curr_x;
        ys_c(add_ct,1) = curr_y;
    end
    
    % Iterate added x,y counter
    add_ct = add_ct + 1;
    
    % Add all indices in the threshold to the consolidated list
    pt_consol(in_thresh) = 1;
    
    % Determine if all points have been consolidated, otherwise find next index
    if sum(pt_consol) == length(pt_consol)
        ALL_PTS_CONSOL = 1;
    else
        curr_ind = find(pt_consol == 0, 1, 'first');
    end
    
end

% Save old points (debugging step)
xs_o = xs;
ys_o = ys;

% Set consolidated points as current
xs = xs_c;
ys = ys_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Order the points using a nearest point algorithm %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select starting point as the bottom-most point
[~,ind_start] = min(ys);

% Initialize the final mask
masked_plane = zeros(size(Xm));
FULL_BND_COMPLETE = 0; % Logical to determine when full boundary has been completed
while FULL_BND_COMPLETE == 0
    % Initialize variables for ordering all points using a nearest neighbor
    % algorithm.
    points_added = zeros(size(ys)); % Points that have been added
    ordered_points = zeros(size(ys)); % Order of points
    ordered_x = zeros(size(ys)); % Ordered x points
    ordered_y = zeros(size(ys)); % Ordered y points
    % Set starting point
    ordered_x(1) = xs(ind_start);
    ordered_y(1) = ys(ind_start);
    ordered_points(1) = ind_start;
    curr_ind = ind_start;
    points_added(ind_start) = 1;
    ct = 2;

    % Iterate through all points until all points have been added
    if iter > 1
        figure(30); hold off; scatter(xs,ys,'ok'); axis([-30 30 -30 30]);
    end
    ALL_POINTS_ADDED = (sum(points_added) == length(ordered_points));
    while ALL_POINTS_ADDED == 0
        % Current x, y position
        curr_x = xs(curr_ind);
        curr_y = ys(curr_ind);

        % Compute distance of all other points
        pt_dist = sqrt((xs-curr_x).^2 + (ys-curr_y).^2);

        % Add large distance for points already added
        pt_dist = pt_dist + 1e3*points_added;

        % Find the next point to add
        [min_dist,next_ind] = min(pt_dist);

        if min_dist > neigh_thresh
            ALL_POINTS_ADDED = 1;
        else
            % Add the next point to the ordered arrays
            ordered_x(ct) = xs(next_ind);
            ordered_y(ct) = ys(next_ind);
            ordered_points(ct) = next_ind;
            points_added(next_ind) = 1;

            % Update next index to current index and update counter
            curr_ind = next_ind;
            ct = ct + 1;

            % Plot the next point
            if iter > 1
                figure(30); subplot(1,2,1); plot(ordered_x(1:ct-1),ordered_y(1:ct-1),'b','LineWidth',1.2); axis([-30 30 -30 30]);
                figure(30); subplot(1,2,2); hold off; scatter(xs,ys,'ok'); axis([-30 30 -30 30]); hold on; 
                plot(xs(points_added == 1),ys(points_added == 1),'or')
                drawnow;
            end

            % Check if all points have been added
            ALL_POINTS_ADDED = (sum(points_added) == length(ordered_points));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Convert the ordered points to the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Need to convert the points (x, y) from mm to pixels
    % Clean up ordered x and y arrays and initialize ordered pixel arrays
    ordered_x = ordered_x(1:ct-1);
    ordered_y = ordered_y(1:ct-1);
    ox_pix = zeros(size(ordered_x));
    oy_pix = zeros(size(ordered_x));

    for q = 1:1:length(ordered_x)
        % Ordered x and y points
        curr_x = ordered_x(q);
        curr_y = ordered_y(q);

        % Get closest index in mask
        [~,x_ind] = min(abs(xp-curr_x));
        [~,y_ind] = min(abs(yp-curr_y));

        % Save closest index
        ox_pix(q) = x_ind;
        oy_pix(q) = y_ind;
    end

    % Get the mask using the polygon identified by the ordered points
    mask_add = double(poly2mask(ox_pix,oy_pix,length(yp),length(xp)));
    
    % Add the mask to the final mask
    masked_plane = masked_plane + mask_add;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Determine if any leftover points need to be revisited %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If the number of points left is greater than 4, find another mask to
    % add to the original mask. Else the mask is complete
    num_pts_left = length(xs) - length(ox_pix);
    if num_pts_left > 4
        % Set the remaining points
        xs = xs(points_added == 0);
        ys = ys(points_added == 0);
        
        % Reset the matrices and starting point
        points_added = zeros(size(ys)); % Points that have been added
        ordered_points = zeros(size(ys)); % Order of points
        ordered_x = zeros(size(ys)); % Ordered x points
        ordered_y = zeros(size(ys)); % Ordered y points
        % Set starting point
        [~,ind_start] = min(ys);
        ordered_x(1) = xs(ind_start);
        ordered_y(1) = ys(ind_start);
        ordered_points(1) = ind_start;
        curr_ind = ind_start;
        points_added(ind_start) = 1;
        ct = 2;
    else
        FULL_BND_COMPLETE = 1;
    end
end



end






