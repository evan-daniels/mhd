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

y2 = zeros(J);
y2 = y + 0.5;

z2 = zeros(K);
z2 = z - 0.5;

rhot = zeros(I,J,K);
rho1 = zeros(J,I,K);
rho2 = zeros(J,I,K);
rho3 = zeros(J,I,K);
rho4 = zeros(J,I,K);

isoval_large = 2.4;
isoval = zeros(2,1);
isoval(1) = 0.5;
isoval(2) = 1.0;

N = input('Enter number of last frame --> ');
M = flag;

for n = 0:N

	count = int2str(n);
	
	s_rho = strcat('rho_files/rho_', count);
	
	if (n < 10)
		count = strcat('0',count);
	end;
	
	S_rho = strcat('full_movie/rho_', count);
	
	rho_vector = load(s_rho);
	
	for i = 1:I
		
		for j = 1:J
		
			for k = 1:K
				rhot(i,j,k) = rho_vector(k + K*((j -1) + (i - 1)*J));
			end;
		end;
	end;
	
	for k = 1:K
		rho1(:,:,k) = rhot(:,:,k)';
	end;
	
	for i = 1:I
		
		for j = 1:J
		
			for k = 1:K
				rho2(j,i,k) = rho1(J+1-j,i,k);
				rho3(j,i,k) = rho1(J+1-j,i,K+1-k);
				rho4(j,i,k) = rho1(j,i,K+1-k);
			end;
		end;
	end;
	
	set(figure,'Colormap',M);
	
	isoval_small = min(min(rho1(:,:,21)));
	
	if (n == 0)
		isoval_small = isoval(1);
	end;
	
	isosurface(x,y,z,rho1,isoval_large);
	isosurface(x,y2,z,rho2,isoval_large);
	isosurface(x,y2,z2,rho3,isoval_large);
	isosurface(x,y,z2,rho4,isoval_large);
	
	axis equal;
	axis([0 1.2 0 1.0 -0.5 0.5]);
	grid on;
	box;
	
	camlight left;
	lighting gouraud;
	alpha(0.2);
	
	isosurface(x,y,z,rho1,isoval_small);
	isosurface(x,y2,z,rho2,isoval_small);
	isosurface(x,y2,z2,rho3,isoval_small);
	isosurface(x,y,z2,rho4,isoval_small);
	lighting gouraud;
	
	print ('-djpeg', '-r300', S_rho);
	
	close;

end;
