% EULER_CLOUD generates the frames for the animation of a shock-cloud
% interaction.
%
% 	CentPack 3d arrays output as single column files (as they are allocated in
% 	memory).  These output is read into matlab as a column vector of length
% 	IxJxK and converted into a matlab 3d array ready to plot.

load('mesh_files/x');
load('mesh_files/y');
load('mesh_files/z');

I = size(x,2);
J = size(y,2);
K = size(z,2);

rhot = zeros(I,J,K);
rho = zeros(J,I,K);

ival = ceil(0.4*K) + 1;
isoval_small = 0.5;
isoval_large = 2.4;

N1 = input('Enter number of first frame --> ');
N2 = input('Enter number of last frame --> ');
M = flag;

for n = N1:N2

	count = int2str(n);
	s_rho = strcat('rho_files/rho_', count);
	
	if (n < 10)
		count = strcat('0',count);
	end;
	
	S_rho = strcat('rho_movie/rho_', count);
	
	rho_vector = load(s_rho);
	
	for i = 1:I
		
		for j = 1:J
		
			for k = 1:K
				rhot(i,j,k) = rho_vector(k + K*((j -1) + (i - 1)*J));
			end;
		end;
	end;
	
	for k = 1:K
		rho(:,:,k) = rhot(:,:,k)';
	end;
	
	set(gcf,'Colormap',M);
  
    isoval_small = 0.5;
	
	isosurface(x,y,z,rho,isoval_large);
	
	axis equal;
	axis([0 1.2 0 0.5 0 0.5]);
	grid on;
	box;
	
	camlight left;
	lighting gouraud;
	alpha(0.2);
	
	isosurface(x,y,z,rho,isoval_small);
	lighting gouraud;
	
	print ('-djpeg', '-r250', S_rho);
	
	close;

end;
