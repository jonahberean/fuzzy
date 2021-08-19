function [x1, x2, x3] = vectortomatrices(Y, N)

%{
    vectortomatrices.m - Takes a 3N^2 dimensional vector, and returns it's representation in terms of three N x N matrices. 
    Inputs:

        N     - coordinate matrix dimension
        Y     - vector to be converted into matrices
                
    Outputs:

        xi    - coordinate matrices

    Other m-files required: none
    Subfunctions:
    MAT-files required: none

    Author: Jonah Berean-Dutcher
    email: jbd@phas.ubc.ca
    Aug 2021; Last revision: 18-Aug-2021 
%}

x = zeros(N,N,N);
count = 1;

% Iterates from N to 1
for d = 1:N

    if(d ~= 1)

        for s = 1:2

            for m = 1:3

                countdiag = 0;
                
                % d = 3, d-1 = 2, (N - (d-1)) = 1
                for i = 1:(N - d + 1)
                    
                    if (s == 1)
                        x(i, d + countdiag,m) = Y(count);
                        count = count + 1; 
                        countdiag = countdiag + 1;
                    else
                        x(d + countdiag,i,m) = Y(count);
                        count = count + 1;
                        countdiag = countdiag + 1;
                    end

                end

            end

        end

    else

        for m = 1:3

            for i = 1:N

                x(i,i,m) = Y(count);
                count = count + 1;

            end
        end

    end


end

x1 = x(:,:,1);

x2 = x(:,:,2);

x3 = x(:,:,3);




