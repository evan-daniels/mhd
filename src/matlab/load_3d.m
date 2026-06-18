% LOAD_3D loads CentPack 3d arrays and converst them into matlab 3d arrays
%
% 	CentPack 3d arrays output as single column files (as they are allocated in
% 	memory).  These output is read into matlab as a column vector of length
% 	IxJxK and converted into a matlab 3d array ready to plot.

function y = load_3d(I,J,K,FILE_NAME)

vart = zeros(I,J,K);
var = zeros(J,I,K);
	
var_vector = load(FILE_NAME);
	
for i = 1:I

    for j = 1:J

        for k = 1:K
            vart(i,j,k) = var_vector(k + K*((j -1) + (i - 1)*J));
        end;
    end;
end;

for k = 1:K
    var(:,:,k) = vart(:,:,k)';
end;

y = var;
    



