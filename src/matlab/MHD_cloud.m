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

pt = zeros(I,J,K);
p = zeros(J,I,K);

N1 = input('Enter number of first frame --> ');
N2 = input('Enter number of last frame --> ');
M = hot(0);

for n = N1:N2

	count = int2str(n);
	
	s_rho = strcat('rho_files/rho_', count);
    s_p = strcat('p_files/p_', count);
	
	if (n < 10)
		count = strcat('0',count);
	end;
	
	S_rho = strcat('rho_movie/rho_', count);
    S_p = strcat('p_movie/p_', count);
	
	rho_vector = load(s_rho);
    p_vector = load(s_p);
	
	for i = 1:I
		
		for j = 1:J
		
			for k = 1:K
				rhot(i,j,k) = rho_vector(k + K*((j -1) + (i - 1)*J));
                pt(i,j,k) = p_vector(k + K*((j -1) + (i - 1)*J));
			end;
		end;
	end;
	
	for k = 1:K
		rho(:,:,k) = rhot(:,:,k)';
        p(:,:,k) = pt(:,:,k)';
	end;
	
	set(gcf,'Colormap',M);
	
	axis equal;
	axis([0 2*pi 0 2*pi 0 2*pi]);
	grid on;
	box;
	
	camlight left;
	lighting gouraud;
	
	slice(x,y,z,rho,0,[],0);
    print ('-djpeg', '-r250', S_rho);
	close;
    
    slice(x,y,z,p,0,[],0);
    print('-djpeg', '-r250', S_p);
	close;

end;
