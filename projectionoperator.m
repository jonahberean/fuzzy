function P = projectionoperator(N, theta)
% projectionoperator.m - generates a projection operator
%
% Syntax: P = projectionoperator(N, theta)
%
% Inputs:
%    N     - matrix dimension, must be odd
%    theta - angle which defines the spherical cap centred at the north
%            pole. The projection operator is meant to project out only
%            those degrees of freedom which are geometrically associated
%            with the region inside the cap
%           
% Outputs:
%    P - the projection operator
%
% Other m-files required:
% Subfunctions:
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% June 2021; Last revision: 14-June-2021
%------------- BEGIN CODE --------------

% for the case of a single matrix (or single scalar field).
antiDiagonalIndex = floor((N/2)*(1 - cos(theta)));
onesVector = ones(1,antiDiagonalIndex);
zerosVector = zeros(1, N-antiDiagonalIndex);
P = diag(horzcat(onesVector, zerosVector));


%% Old, from the attempt to do all three fluctuation matrices.
% % When we construct K, we arrange the x_j fluctuation matrix elements into
% % a column vector, previously named Y. We need to identify which of the
% % components of Y, the y_i, correspond to matrix elements [x_c]_{a,b} for
% % which (a + b) < N(1 - cos(theta)). With those indices i identified, we
% % can construct the appropriate projection operator within the 3N^2
% % dimensional space of the vector Y. 
% 
% % initialize a matrix, P, of zeros that we will fill with the non-zero
% % entries of the projection operator.
% P = zeros(3*N^2, 3*N^2);
% 
% for a = 1:N
%     
%     for b = 1:N
%         
%         if ((a + b) < N*(1-cos(theta)))
%             
% 
%             % iterate over the three fluctuation matrices
%             for c = 1:3
% 
%                 d = b - a;            
%                 I = + abs(sign(d)) * 3 * (N + 2 * (N * (abs(d) - 1) - ((abs(d)-1)*(abs(d)))/2))...
%                     + abs(min(sign(d),0)) * (3 * (N - abs(d))) ...
%                     + (c-1) * (N - abs(d)) ...
%                     + (a + min(0,d)); 
% 
%                 P(I,I) = 1;
% 
%             end
%         
%         end        
%         
%     end
%     
% end