function [u,v,p] = showData (file)
%
%  Just write inside the Matlab-shell: 
%     showData('./../<filename.out>')
%
fid = fopen(file);
xlength = fread(fid,1,'float');
ylength = fread(fid,1,'float');
imax = fread(fid,1,'float');
jmax = fread(fid,1,'float');
u = ones(imax,jmax);
v = u;
p = u;

u    = fread(fid,[imax,jmax],'float');
v    = fread(fid,[imax,jmax],'float');
p    = fread(fid,[imax,jmax],'float');
temp    = fread(fid,[imax,jmax],'float');
zeta    = fread(fid,[imax,jmax],'float');
psi    = fread(fid,[imax,jmax],'float');

fclose(fid);

% Open one window 
%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[100 100 700 700])

% Draw the level sets of u,v,p,psi,... inside one window 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position',[.05 .55 .4 .4])  % coordinates inside the window 
contour(u',40);
title('Isolines of u-velocity')

axes('position',[.05 .05 .4 .4])
contour(v',40);
title('Isolines of v-velocity')

%axes('position',[.55 .55 .4 .4])
%contour(p',40);
%title('Isolines of pressure')

axes('position',[.55 .55 .4 .4])
contour(psi',40);
axis([0 imax 0 jmax]);
title('stream function Psi')

axes('position',[.55 .05 .4 .4])
quiver(u',v',5);
%axis([0 50 0 50]);
axis([0 imax 0 jmax]);
title('Velocity vector field')

return
