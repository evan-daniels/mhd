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

ut = zeros(I,J,K);
u = zeros(J,I,K);

isoval_large = 0.5;
isoval = zeros(2,1);
isoval(1) = 0.1;
isoval(2) = 0.5;

N1 = input('Enter number of first frame --> ');
N2 = input('Enter number of last frame --> ');
M = flag;

for n = N1:N2

	count = int2str(n);
	
	s_u = strcat('u_files/u_', count);
	
	if (n < 10)
		count = strcat('0',count);
	end;
	
	S_u = strcat('u_movie/u_', count);
	
	u_vector = load(s_u);
	
	for i = 1:I
		
		for j = 1:J
		
			for k = 1:K
				ut(i,j,k) = u_vector(k + K*((j -1) + (i - 1)*J));
			end;
		end;
	end;
	
	for k = 1:K
		u(:,:,k) = ut(:,:,k)';
	end;
	
	set(gcf,'Colormap',M);
	
	isoval_small = 0.1;
	
	isosurface(x,y,z,u,isoval_small);
	
	axis equal;
	axis([0 1.0 0 1.0 0 1.0]);
	grid on;
	box;
	
	camlight left;
	lighting gouraud;
	alpha(0.2);
	
	isosurface(x,y,z,u,isoval_large);
	lighting gouraud;
	
	print ('-djpeg', '-r250', S_u);
	
	close;

end;
