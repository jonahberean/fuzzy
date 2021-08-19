function mu = findmu(N, xi, xi2)
% findmu.m - identifies the critical values of mu at given N, and for given matrices xi and xi2
%
% Inputs:
%    N    - coordinate matrix dimension
%    xi   - projected height operator into some subspace of the full 3N^2 dimensional matrix theory
%    xi2  - same projection, but for the square of the height operator
%
% Outputs:
%    mu   - array of critical mu values
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% May 2021; Last revision: 29-Jul-2021
%------------- BEGIN CODE --------------

% round the reduced xi and xi2 matrices until they are symmetric
roundToDigit = 18;
tf = false;
while (~tf)
    xi = round(xi, roundToDigit);
    xi2 = round(xi2, roundToDigit);
    tf = ishermitian(xi) & ishermitian(xi2);
    roundToDigit = roundToDigit - 1;
end

% define array of values of mu to search over
eigvals = eig(xi);
muStep = var(eigvals(1:3));
fprintf('The eigenvalues of xi:\n')
eigvals
muSearchArray = 1:muStep:2*(N-2);
fprintf('Will search over %d values of mu in the range[%d to %d]\n', length(muSearchArray), muSearchArray(1), muSearchArray(end))
tolerance = muStep / 100;

% define arrays to hold the lowest two eigenvalues, their mean, and
% their variance
lowestEigval    = zeros(1, length(muSearchArray));

% iterate over mu
fprintf('Iterating over mu.\n')
for i = 1:length(muSearchArray)

%     if (mod(i,100) == 0)
% 
%         % print status to console
%         statement = 'Searching over mu: %d out of %d\n';
%         fprintf(statement, i, length(muSearchArray)) 
% 
%     end

    % diagonalize the search matrix at this mu value
    eigvals = eig(xi2 - muSearchArray(i)*xi);

    % retrieve the lowest two eigenvalues
    lowestEigval(i) = eigvals(1);

end

% array of the simult eigvec slopes
bArray = mod(N,2):2:(N-2);

overlapMu = zeros(1, length(bArray)*2);

% We should be able to retrieve all the array elements where they overlap
for i = 1:length(bArray)
    
    line = bArray(i)^2 - muSearchArray.*bArray(i);
    
    % identify points where the absolute difference between this line and
    % the lowestEigval array falls below tolerance
    diff = abs(line - lowestEigval);
    muRange = muSearchArray(diff < tolerance);
    
    if (~isempty(muRange))
        
        overlapMu(2*i - 1) = muRange(1); 
        overlapMu(2*i)     = muRange(end);
        
    end
     
end


% mu = 1;
mu = overlapMu;