function Z = heightoperator(N)
% heightoperator.m - Generates the height operator for 3N^2 dim. system
%
% Syntax:  Z = heightoperator(N)
%
% Inputs:
%    N - coordinate matrix dimension
%
% Outputs:
%    Z - height operator for 3N^2 dim. system
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 16-March-2021
%------------- BEGIN CODE --------------

% % A symbolic approach to make sense of what we're doing.

% % returns the three SU(2) generators in the N-dim. irreducible rep. 
% [L1, L2, L3, ~, ~] = su2generators(N);

% % Three NxN fluctuation matrices containing symbolic variables, 
% % real valued.
% syms Rx1 [N N];
% syms Rx2 [N N];
% syms Rx3 [N N];

% % Substituting in the real valued matrix elements
% x1 = 1/2*(Rx1 + Rx1.') - 1/(2*sqrt(-1))*(Rx1 - Rx1.');
% x2 = 1/2*(Rx2 + Rx2.') - 1/(2*sqrt(-1))*(Rx2 - Rx2.');
% x3 = 1/2*(Rx3 + Rx3.') - 1/(2*sqrt(-1))*(Rx3 - Rx3.');

% % Returns a column vector containing the matrix elements of the three real-valued fluctuation matrices
% Y = matricestovector(Rx1, Rx2, Rx3, N);

% % Acts on the complex values fluctuation matrices with the Z operator, where the x1, x2, and x3 in these expressions (and hence the resulting Zx1, Zx2, Zx3 also) are expressions in terms of the real values matrices.
% Zx1 = (1/2)*(x1*L3 + L3*x1);
% Zx2 = (1/2)*(x2*L3 + L3*x2);
% Zx3 = (1/2)*(x3*L3 + L3*x3);

% % Now we convert back from the complex space to the real space and get the effect of the Z action upon the real-valued matrices.
% RZx1 = 1/2*(Zx1+Zx1.')+1/(2*sqrt(-1))*(Zx1-Zx1.');
% RZx2 = 1/2*(Zx2+Zx2.')+1/(2*sqrt(-1))*(Zx2-Zx2.');
% RZx3 = 1/2*(Zx3+Zx3.')+1/(2*sqrt(-1))*(Zx3-Zx3.');

% % Generates the column vector with entries that are expressions of real-valued matrix elements, expressions that have been transformed by acting with Z on the complex valued space.
% ZY = matricestovector(RZx1, RZx2, RZx3, N);

% % Construct the matrix that carries out the Y -> ZY transformation by
% % taking partial derivatives
% Z = zeros(3*N^2, 3*N^2);

% for i = 1:3*N^2

%     fprintf("i = %d out of %d\n", i, 3*N^2)

%     for j = 1:3*N^2

%         Z(j,i) = diff(ZY(j), Y(i), 1);

%     end

% end

% Znew = Z;

%% Older code without proper rescaling to make the sphere unit-valued:    

% initializing a vector of zeros which will become the main diagonal of Z
Z = zeros(3*N^2,1);

counter = 1;

% looping over the values of |l| = 0, 1, ... N-1, which index the diagonals of the
% fluctation matrices
for l = 0:(N-1)
    
    % loop that corresponds to iterating over the three fluctuation matrices
    for j = 0:2
        
        % In the lth diagonal, there are N-|l| elements. Here we loop over the 
        % elements of this diagonal, from 1 to N-|l|
        
        for k = 0:(N-l-1)
            
            Z(counter) = (N-l-1-2*k);
            counter = counter + 1;
            
        end
        % if it's not the main diagonal, then we repeat the same elements for the (-l)
        % diagonal
        if ((l ~= 0) && (j == 2))
            
            
            Z(counter:counter + (3*(N-l)-1)) = Z(counter - 3*(N-l):counter-1);
            counter = counter + 3*(N-l);
            
        end
        
    end
end

% Normalizing to make the sphere's radius unit valued.
nu = 2 / sqrt(N^2 - 1);

% As per the definition, (1/2)(x_i L_3 + L_3 x_i)
Z = (1/2)*nu*diag(Z);