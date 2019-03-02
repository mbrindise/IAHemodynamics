% Number of lines by data file:
% ANY = 111 - NAME = Tracks_Run1_Full.dat - # Lines = 72,695,624
% ANY = 007 - NAME = Tracks_Run1_Full.dat - # Lines = 18,899,827
% ANY = 007 - Name = Tracks_Run1_1to1800.dat - # Lines = 4,616,039
% ANY = 007 - Name = Tracks_Run1_1750to2050.dat - # Lines = 618,452
% ANY = 007 - Name = Tracks_Run1_2000toEnd.dat - # Lines = 11,304,983
% ANY = 007 - Name = Tracks_Run2_1to4000.dat - # Lines = 20,334,513
% ANY = 007 - Name = Tracks_Run2p1.dat - # Lines = 65,494,138
% ANY = 007 - Name = Tracks_Run3p1.dat - # Lines = 76,032,478
% Or 007 - Name = Track_Orientation_ANY007_Horiz_Run02 - # Lines = 43,004,036
% Or 007 - Name = Track_Orientation_ANY007_Vert_Run02 - # Lines = 44,998,481
% Or 111 - Name = Track_Orientation_Horiz_Run04p2 - # Lines = 20,546,848
% Or 111 - Name = Track_Orientation_Vert_Run04p2 - # Lines = 22,017,148

num_lines = 43004036;
NUM_LINES_UNKNOWN = 0;
start_time = 1;
t_start = 1;
t_end = 5766;
% File path of the main Shake the Box file
any_num = '007';
cpu_sys = 'pc'; % Can be server of PC

if strcmp(cpu_sys,'pc')
    main_file_path = ['Z:\Projects\Cerebral_Aneurysm\ANY\',any_num,'\STB\'];
    fslash = '\';
else
    main_file_path = ['~/Projects/Cerebral_Aneurysm/ANY/',any_num,'/STB/'];
    fslash = '/';
end

% Name of Shake the Box file which contains all track across all time
main_file_name = 'Tracks_Orientation_ANY007_Horiz_Run02.dat';

% Load the big file which contains all data points. This file is often far
% too big to open in full (requires too much RAM so only a part can be
% opened at a time)
fid = fopen([main_file_path,main_file_name],'r');
if NUM_LINES_UNKNOWN
    % Open the first column only
    A = textscan(fid,'%f %*[^\n]','HeaderLines',2);
    num_lines = length(A{1,1});
else
    % Load in all point values
    xa = nan(num_lines,1);
    ya = nan(num_lines,1);
    za = nan(num_lines,1);
    ua = nan(num_lines,1);
    va = nan(num_lines,1);
    wa = nan(num_lines,1);
    time_step = nan(num_lines,1);
    track_id = nan(num_lines,1);
    bad_time_step = 0;
    % Print for information
    fprintf('\r\nLoading all track data... ')
    prev_done = 0;
    ct = 1;
    for q = 1:1:num_lines
        % STB text files go: [x, y, z, u, v, w, |v|, I, time step, track ID, ax, ay, az, |a|]
        A = textscan(fid,'%f %f %f %f %f %f %*f %*f %f %f %*f %*f %*f %*f',1,'HeaderLines',2*(q==1),'Delimiter',' ');
        %A = textscan(fid,'%f',1,'HeaderLines',1+q,'Delimiter',' ');
        if ~isempty(A{1,1})
            xa(q) = A{1,1};
            ya(q) = A{1,2};
            za(q) = A{1,3};
            ua(q) = A{1,4};
            va(q) = A{1,5};
            wa(q) = A{1,6};
            time_step(q) = A{1,7};
            track_id(q) = A{1,8};
        else
            bad_time_step(ct) = q;
            ct = ct + 1;
            keyboard;
        end
        
        % Print the percent done
        perc_comp = floor(100*q/num_lines);
        if mod(perc_comp,5) == 0 && perc_comp ~= prev_done
            fprintf(' %i%%',perc_comp)
            prev_done = perc_comp;
        end
    end
end
fclose(fid);
fprintf(' complete\r\n')

% Seperate points by time step and save time steps
if ~NUM_LINES_UNKNOWN
    fp_out_dir = [main_file_path,'raw',fslash];
    fp_out_name = ['any',any_num,'_',any_fp_orient,'_stb_v1_'];
    unq_time = unique(time_step) + start_time - 1;
    if isempty(t_start), t_start = unq_time(1); end
    if isempty(t_end), t_end = unq_time(end); end
    fprintf('\r\n')
    for curr_time = t_start:1:t_end
        curr_time_ind = curr_time - start_time + 1;
        % Get all points in the current time
        time_pts = find(time_step == curr_time_ind);
        x = xa(time_pts);
        y = ya(time_pts);
        z = za(time_pts);
        u = ua(time_pts);
        v = va(time_pts);
        w = wa(time_pts);
        trackID = track_id(time_pts);

        save([fp_out_dir,fp_out_name,num2str(curr_time,'%06i'),'.mat'],'x','y','z','u','v','w','trackID');
        
        % Save dat file
        if curr_time < 5
            % Initialize output file
            fid = fopen([fp_out_dir,fp_out_name,num2str(curr_time,'%06i'),'.dat'],'w');
            fprintf(fid,['TITLE = "',fp_out_name,num2str(curr_time,'%05i'),'"']);
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
        fprintf('Completed: %i \n',curr_time)
    end
end