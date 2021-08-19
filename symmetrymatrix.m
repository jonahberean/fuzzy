function B = symmetrymatrix(N)
% symmetrymatrix - Generates matrix representation of so(3) element
%
% Syntax:  B = symmetrymatrix(N)
%
% Inputs:
%    N - coordinate matrix dimension
%
% Outputs:
%    B - matrix representation of so(3) operator acting on 3N^2 dim system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 16-March-2021
%------------- BEGIN CODE --------------


% B0 has a separate structure from B(N-1), we construct it first
B0 = [zeros(N) -eye(N)  zeros(N); ...
      eye(N)   zeros(N) zeros(N); ...
      zeros(N) zeros(N) zeros(N)];
  
% !!! Supposing I've got it wrong, switching sign:
% B0 = [zeros(N) eye(N)  zeros(N); ...
%       -eye(N)   zeros(N) zeros(N); ...
%       zeros(N) zeros(N) zeros(N)];
  
B = blkdiag(B0);

% Looping from 0 to N-1, constructing the blocks of B
for l = 1:N-1
    
    Bl = [zeros(N-l) -eye(N-l)  zeros(N-l) -l*eye(N-l) zeros(N-l)  zeros(N-l); ... 
          eye(N-l)   zeros(N-l) zeros(N-l) zeros(N-l)  -l*eye(N-l) zeros(N-l); ...
          zeros(N-l) zeros(N-l) zeros(N-l) zeros(N-l)  zeros(N-l)  -l*eye(N-l); ...
          l*eye(N-l) zeros(N-l) zeros(N-l) zeros(N-l)  -eye(N-l) zeros(N-l); ...
          zeros(N-l) l*eye(N-l) zeros(N-l) eye(N-l)    zeros(N-l)  zeros(N-l); ...
          zeros(N-l) zeros(N-l) l*eye(N-l) zeros(N-l)  zeros(N-l)  zeros(N-l)];
    
    % Adding the new block to the block diagonal B as we go 
    B = blkdiag(B, Bl);
end

