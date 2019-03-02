% Number of lines by data file:
% ANY = 111 - Time = 1.623000 - 0.0015 - 2.434500
% ANY = 007 - Time = 1.161000 - 0.0015 - 1.74000
% ANY = 20160926 - Time = 1.278 - 0.002 - 1.286
% ANY = 20170825 - Time = 1.384 - 0.002 - 1.392
% ANY = 20170906 - Time = 1.132 - 0.002 - 1.140
% ANY = 20170911 - Time = 1.534 - 0.002 - 1.542
% ANY = 20171020 - Time = 1.490 - 0.002 - 1.498

t_start = 1.534;%1.623; %1.161000;%1.623000;
dt = .00200;
t_end = 1.542;%2.4345;% 1.740000;%2.434500;
% File path of the main Shake the Box file
any_num = '20170911';
main_file_path = ['Z:\Projects\Cerebral_Aneurysm\ANY\',any_num,'\CFD\fluent_transfer\'];
fp_out_dir = ['Z:\Projects\Cerebral_Aneurysm\ANY\',any_num,'\CFD\raw\'];
fp_out_name = ['any',any_num,'_cfd_'];
    
% Name of Shake the Box file which contains all track across all time
base_file_name = '20170911CFD_time';% 'UCSF_ANY111_150_Transient-';%';%'UCSF_ANY111_150_Transient-';%
file_ext = '.dat';% 'Test_Tracks.dat';

% Load the big file which contains all data points. This file is often far
% too big to open in full (requires too much RAM so only a part can be
% opened at a time)
ct = 1;
t = 0;
num_iter = length(t_start:dt:t_end);
iter = 1;
prev_done = 0;
fprintf('\nSaving CFD raw files...')
for q = t_start:dt:t_end
    A = importdata([main_file_path,base_file_name,num2str(q,'%.3f'),file_ext],',');
    A = A.data;
    % Obtain the info
    if ~isempty(A(1,1))
        x = 1e3*A(:,2);
        y = 1e3*A(:,4);
        z = 1e3*A(:,3);
        u = A(:,6);
        v = A(:,8);
        w = A(:,7);
        %wssX = A(:,11);
        %wssY = A(:,12);
        %wssZ = A(:,13);
    else
        bad_time_step(ct) = q;
        ct = ct + 1;
        keyboard;
    end
    
    % Save mat file
    save([fp_out_dir,fp_out_name,num2str(iter,'%05i'),'.mat'],'x','y','z','u','v','w','t');
    %save([fp_out_dir,fp_out_name,'wss_',num2str(iter,'%05i'),'.mat'],'wssX','wssY','wssZ');
    
    % Save dat file - only for first 10 because of time/file size of
    % original CFD dat files
    if iter <= 1
        % Initialize output file
        fid = fopen([fp_out_dir,fp_out_name,num2str(iter,'%05i'),'.dat'],'w');
        fprintf(fid,['TITLE = "',fp_out_name,num2str(t,'%05i'),'"']);
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
    % Print the percent done
    perc_comp = floor(100*iter/num_iter);
    if mod(perc_comp,5) == 0 && perc_comp ~= prev_done
        fprintf(' %i%%',perc_comp)
        prev_done = perc_comp;
    end
    % Increment the iteration
    iter = iter + 1;
    % Increment the time
    t = t + dt;
end

fprintf(' complete\r\n')


