function [uf,vf,wf,eval] = UOD3(u,v,w,mask,winsize,tol,minvel,mag)
%%% Universal outlier detection code in 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: u, v, w : velocity components, each entered as 3D arrays
%         mask : should be same size as velocity field and indicates where
%         the valid velocity points exist in the grid. (1 for valid, 0 for
%         invalid)
%         winsize: Window size (will be a cube, only one number needs to be
%                  input. Takes first number if more than one is input)
%         tol: Tolerance of the median deviation (typical threshold: 2-3)
%         minvel: Minimum number of good velocity points for UOD to be
%         valid
%         mag: magnification to convert eps to grid. This should be the dx
%         of the velocity field
%
% OUTPUTS: uf, vf, wf : velocity arrays that have been validated using UOD
%          eval: logical array of same size as velocity fields that 
%                indicates if the value at that point was replaced
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate inputs
% Get size of velocity field
[nR,nC,nZ] = size(u);
ws = winsize(1,1); % Obtain window size
chalfws = ceil(ws/2); % Half window size
fhalfws = floor(ws/2); % Half window size

%% Set constants
e = 0.1*mag; % Minimum variance (in m)

%% Initialize output
eval = nan(size(u));
uf = nan(size(u));
vf = nan(size(v));
wf = nan(size(w));


%% Perform UOD for all points (except boundaries)
for k = chalfws:1:nZ-chalfws
    for j = chalfws:1:nC-chalfws
        for i = chalfws:1:nR-chalfws
            % Get current velocity at point
            uc = u(i,j,k);
            vc = v(i,j,k);
            wc = w(i,j,k);
            % Get current velocity blocks
            ublock = u(i-fhalfws:i+fhalfws,j-fhalfws:j+fhalfws,k-fhalfws:k+fhalfws);
            vblock = v(i-fhalfws:i+fhalfws,j-fhalfws:j+fhalfws,k-fhalfws:k+fhalfws);
            wblock = w(i-fhalfws:i+fhalfws,j-fhalfws:j+fhalfws,k-fhalfws:k+fhalfws);
            
            % Set the current velocity of the block as NaN so it does not
            % bias the median calculation
            ublock(chalfws,chalfws,chalfws) = NaN;
            vblock(chalfws,chalfws,chalfws) = NaN;
            wblock(chalfws,chalfws,chalfws) = NaN;
            % Set the masked portions to NaN
            mblock = mask(i-fhalfws:i+fhalfws,j-fhalfws:j+fhalfws,k-fhalfws:k+fhalfws);
            ublock(mblock == 0) = NaN;
            vblock(mblock == 0) = NaN;
            wblock(mblock == 0) = NaN;
            gdpts = ~isnan(ublock);
            numgd = sum(gdpts(:));
            
            % Determine the number of good points
            if numgd >= minvel
            
                % Compute median of block
                umed = nanmedian(ublock(:));
                vmed = nanmedian(vblock(:));
                wmed = nanmedian(wblock(:));

                % Compute residual and redisual median
                ures = abs(ublock - umed);
                vres = abs(vblock - vmed);
                wres = abs(wblock - wmed);
                rumed = nanmedian(ures(:));
                rvmed = nanmedian(vres(:));
                rwmed = nanmedian(wres(:));

                % Compute residual of current point
                ru = abs(uc - umed)/(rumed + e);
                rv = abs(vc - vmed)/(rvmed + e);
                rw = abs(wc - wmed)/(rwmed + e);
                restot = sqrt(ru.^2 + rv.^2 + rw.^2); % Total residual

                % Determine if residual is greater than the tolerance, if so,
                % set the current value as the median, otherwise, set it to the
                % input - the value passes
                if restot > tol
                    uf(i,j,k) = umed;
                    vf(i,j,k) = vmed;
                    wf(i,j,k) = wmed;
                    eval(i,j,k) = 1;
                else
                    uf(i,j,k) = uc;
                    vf(i,j,k) = vc;
                    wf(i,j,k) = wc;
                end
            end
        end
    end
end



end

