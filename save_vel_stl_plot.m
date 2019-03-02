function [] = save_vel_stl_plot(cadat)

% Use the first time step
ts = cadat.options.TS;

% Load the velocity file
fp = load([cadat.DIR.velfiles,cadat.FILE.basename,num2str(ts,'%05i'),'.mat']);
% Get the velocity field
x = fp.x; y = fp.y; z = fp.z;
u = fp.u; v = fp.v; w = fp.w;
% Convert the current field to a mesh grid
% Get dx, dy, dz
unqV = unique(x);
dx = unqV(2)-unqV(1);
unqV = unique(y);
dy = unqV(2)-unqV(1);
unqV = unique(z);
dz = unqV(2)-unqV(1);
% Create lines
xl = min(x)-2*dx:dx:max(x)+2*dx;
yl = min(y)-2*dy:dy:max(y)+2*dy;
zl = min(z)-2*dz:dz:max(z)+2*dz;
% Create mesh
[X,Y,Z] = meshgrid(xl,yl,zl);
U = zeros(size(X)); V = zeros(size(X)); W = zeros(size(X));
% Place velocities on the mesh grid
for zz = 1:1:length(x)
    % Current velocity point
    cx = x(zz); cy = y(zz); cz = z(zz);
    cu = u(zz); cv = v(zz); cw = w(zz);
    % Find the minimum distance
    [~,r_ind] = min(abs(yl-cy));
    [~,c_ind] = min(abs(xl-cx));
    [~,z_ind] = min(abs(zl-cz));
    % Place the velocity
    U(r_ind,c_ind,z_ind) = cu;
    V(r_ind,c_ind,z_ind) = cv;
    W(r_ind,c_ind,z_ind) = cw; 
end
% Close figure(123) if it is still open
figure(123); close(gcf);
% Draw surface
vmag = sqrt(U.^2 + V.^2 + W.^2);
vv = double(vmag > 0);
figure(123); hold off;
p = patch(isosurface(X,Y,Z,vv,0));
p.FaceColor = [0.7 0.7 0.7];
p.EdgeColor = 'none';
p.FaceAlpha = 0.2;


% Save stl of the geometry for plotting
faces = p.Faces;
verts = p.Vertices;
stl_out_name = [cadat.DIR.stlfiles,cadat.FILE.anyname,'velocity_plot','.stl'];
stlwrite(stl_out_name,faces,verts)


end

