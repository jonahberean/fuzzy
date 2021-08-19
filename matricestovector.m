function Y = matricestovector(x1, x2, x3, N)

%region - doc
%{
matricestovector.m - Takes three symbolic N x N matrices, and returns the
entries in a column vector

Inputs:
    N     - coordinate matrix dimension
           
Outputs:
    Y - column vector


Other m-files required: 
Subfunctions:
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
Aug 2021; Last revision: 18-Aug-2021
%}
%endregion - doc

% Create an empty vector and iteratively append diagonals from the input matrices to it
Y = sym(zeros(1,3*N^2));

% Append the main diagonals to Y
Y(1:N)       = diag(x1);
Y(N+1:2*N)   = diag(x2);
Y(2*N+1:3*N) = diag(x3);

% index of the first element after the main diagonals
start = 3*N+1;

% placeholder for counting indices
step = 0;

% Looping from 1 to N-1, the number of diagonals of a given sign i.e. 
% upper or lower
for i = 1:N
    start = start + step;
%     disp(start)
    Y(start           : start+(N-1-i))         = diag(x1, i);
    Y(start+1*(N-i)   : start+2*(N-i)-1)       = diag(x2, i);
    Y(start+2*(N-i)   : start+3*(N-i)-1)       = diag(x3, i);
    Y(start+3*(N-i)   : start+4*(N-i)-1)       = diag(x1, -i);
    Y(start+4*(N-i)   : start+5*(N-i)-1)       = diag(x2, -i);
    Y(start+5*(N-i)   : start+6*(N-i)-1)       = diag(x3, -i);

    step = 6*(N-i);
end

