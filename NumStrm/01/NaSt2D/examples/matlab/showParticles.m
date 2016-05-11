function  [p]=showParticles(name);
%
% Reads and visualizes all particles from "<name>.tra" or "<name>.str" files.
%
%  Just write inside the Matlab-shell:
%     showParticles('./../<filename.str>')
%   or
%     showParticles('./../<filename.tra>')
%
% Reads the positions (x,y) of the particles from file 'name'.
% The positions are scaled and rounded to fit into a 500x500 sparse matrix p.
% The matrix p is returned.
%

colors = 'brgymck';

[fid,msg]=fopen(name,'r') ;

if (fid<=0) 
   msg 
    return ;
  end

  % Read the stuff 
  [imax, c]=fread(fid,1,'int') ;
  [jmax, c]=fread(fid,1,'int') ;
  [delx, c]=fread(fid,1,'double') ;
  [dely, c]=fread(fid,1,'double') ;
  [lines, c]=fread(fid,1,'int') ;

  xl(1)=imax*delx;
  xl(2)=jmax*dely;

  [n ,c]=fread(fid,1,'float');
  while n>1,
    for i=0:(lines-1),
      x=zeros(n,2) ;
      % Scale in each direction
      xm=min(xl(1),xl(2)) ;
      nx=round(500*xl(1)/xm) ;
      ny=round(500*xl(2)/xm) ;
      % Read the positions
      [x,c]=fread(fid,[2,n],'float') ;
      x=x' ;
      % round and scale
      x(:,1)=max(floor((nx*x(:,1)/xl(1))+1), 1);
      x(:,2)=max(floor((ny*x(:,2)/xl(2))+1), 1);
      n = length(x(:,1));
      part=sparse(x(:,1),501-x(:,2) , ones(n,1) ,nx+1 , ny+1);
      spy(part', sprintf('-.%s', colors(mod(i, 7) + 1)), 1);
      title('PRESS SPACE');
      hold on;
      [n ,c]=fread(fid,1,'float');
    end
    pause;
    hold off;
  end  
  title END;
  fclose(fid) ;
						  
