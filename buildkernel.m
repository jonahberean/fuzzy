function nullEigenvectors = buildkernel(N)
% buildkernel.m - diagonalizes the coupling matrix and returns the
% eigenvectors spanning the kernel, saves the result to file
%
%
% Syntax: V = buildkernel(K)
%
% Inputs:
%    N - coordinate matrix dimension
%
% Outputs:
%    V - matrix whose columns are null eigenvectors of K
%
% Other m-files required:
% Subfunctions:
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% April 2021; Last revision: 27-April-2021
%------------- BEGIN CODE --------------

filename = ...
        [fileparts(pwd), '/pythonK/direct_', num2str(N), '.mat'];
M = load(filename); 
K = M.K;

if(~issymmetric(K))
    disp('Warning, K not symmetric')
else
    disp('K is symmetric')
end

% diagonalizes K
[V, D] = eig(K);

% Retrieves and returns just the null eigenvectors
nullEigenvectors = V(:, abs(diag(D)) < 1e-4);

% Saves the nullEigenvectors to file
save([fileparts(pwd), '/savedkernel/', num2str(N),'.mat'],...
    'nullEigenvectors');
