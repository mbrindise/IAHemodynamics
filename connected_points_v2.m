function [point_segment] = connected_points_v2(x,y,pt_dist_thresh)

num_pts = length(x);
dist_thresh = pt_dist_thresh; % Distance threshold

% Initialize variables
point_checked = zeros(size(x)); % Create points connected vector which determines which points have been completed
point_segment = zeros(size(x)); % Create array to store what segment each point belongs to

% Iterate through all points until all points have been connected to a
% segment
iter = 0;
curr_ind = 1;
curr_seg = 1;
point_segment(curr_ind) = curr_seg;
point_checked(curr_ind) = 1;
ALL_POINTS_COMPLETE = 0;
while ALL_POINTS_COMPLETE == 0
    SEGMENT_COMPLETE = 0;
    inds2chk = [];
    while SEGMENT_COMPLETE == 0
        curr_pt_x = x(curr_ind);
        curr_pt_y = y(curr_ind);

        % Check if points are within the distance threshold
        dist_pts = sqrt((x - curr_pt_x).^2 + (y - curr_pt_y).^2);    
        dist_pass = dist_pts <= dist_thresh; % Check if the point is within the distance threshold
        new_pt = dist_pts ~= 0; % Ensure the point is not itself
        new_pt2 = point_segment == 0; % Ensure the point has not been previously assigned
        new_log = dist_pass + new_pt + new_pt2;
        new_inds = find(new_log == 3);
        
        % Add the new indices to the current segment
        point_segment(new_inds) = curr_seg;
        
        % Add the new unchecked indices as new points to check for additional
        % neighbors within the segment
        for zz = 1:1:length(new_inds)
            c_new_ind = new_inds(zz);
            % Check if the current new index has been checked, if not, add
            % it to the list
            if point_checked(c_new_ind) == 0
                inds2chk = [inds2chk,c_new_ind];
            end
        end

        % Update the current index and remove it from list to check. If the
        % list of indices to check is empty then the segment is complete.
        if ~isempty(inds2chk)
            curr_ind = inds2chk(1); 
            point_checked(curr_ind) = 1; % Mark the new current index as checked
            inds2chk = inds2chk(2:end); % Remove it from indices to check
        else
            SEGMENT_COMPLETE = 1;
            curr_seg = curr_seg + 1;
        end
        iter = iter + 1;
    end
    
    if sum(point_segment == 0) == 0
        ALL_POINTS_COMPLETE = 1;
    else
        curr_ind = find(point_segment == 0, 1, 'first');
        point_segment(curr_ind) = curr_seg;
    end
end

end




