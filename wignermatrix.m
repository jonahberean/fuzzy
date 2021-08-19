function W = wignermatrix(N, beta)
% wignermatrix.m -
%
% Syntax: wignermatrix(N)
%
% Inputs:
%    N     - coordinate matrix dimension
%           
% Outputs:
%
% Other m-files required: wignermatrix.m
% Subfunctions:
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% Aug 2021; Last revision: 10-Aug-2021

j = (1/2)*(N-1);

% vector of m values 
mValues = j:-1:-j;
% instantiate the matrix itself
W = zeros(N);

% loop over the row index
for i = 1:N
    
    % compute m (m in Sakurai's notation)
    m = mValues(i);
    
    % loop over the column index
    for ii = 1:N
        
        % compute n (m' in Sakurai's notation)
        n = mValues(ii);
        
        % instantiate the matrix element
        element = 0;
        
        % loop over k from 0 to j + m
        for k = 0:2*j
            
            sig = k - n + m;
            n1 = factorial(j + n);
            n2 = factorial(j - n);
            n3 = factorial(j + m);
            n4 = factorial(j - m);

            d1 = j + n - k;
            d2 = k;
            d3 = j - k - m;
            d4 = k - n + m;

            denominatorTest = d1 < 0 || d2 < 0 || d3 < 0 || d4 < 0;

            % any of the factorial arguments appearing in the denominator are
            % negative, then we do not compute the term in the sum for this value
            % of k
            if(denominatorTest)

                result = 0;

            % else we do compute it
            else

                % compute factorials of the four denominator factors
                d1 = factorial(d1); 
                d2 = factorial(d2); 
                d3 = factorial(d3); 
                d4 = factorial(d4);

                % compute numerator and denominator
                num = sqrt(n1 * n2 * n3 * n4);
                den = d1 * d2 * d3 * d4;

                % exponents on the cosine and sine
                cpow = 2*j - 2*k + n - m;
                spow = 2*k - n + m;

                result = (-1)^sig * (num / den) * cos(beta/2)^cpow * ...
                    sin(beta/2)^spow;

            end  
            
            % add the result to the matrix element
            element = element + result;
            
        end
        
        % set the value of the matrix element in the matrix
        W(i,ii) = element;
        
    end
    
end
    