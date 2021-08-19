function Z = getHeightOperator(N)  

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