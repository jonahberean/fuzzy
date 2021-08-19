function [T1,T2,T3,Tminus,Tplus] = su2generators(N)
% FUNCTION_NAME - returns the SU(2) generators in the N-dim. irrep. 
%
% Syntax:  [T1,T2,T3,Tminus,Tplus] = su2generators(N)
%
% Inputs:
%    N - coordinate matrix dimension
%
% Outputs: !!!
%    output1 - Description
%    output2 - Description
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% March 2021; Last revision: 16-March-2021
%------------- BEGIN CODE --------------
%% su2generators

    j = (N-1)/2;

    mvalues = -j:1:j;

    % Tminus is zeros, except for the first diagonal below the main diagonal
    Tminus = zeros(length(mvalues));
    % Tminus is zeros, except for the first diagonal above the main diagonal
    Tplus  = zeros(length(mvalues));

    % Generating Tplus
    % iterating over the elements of the (-1th) diagonal
    for i=2:(length(mvalues))
        p = i-1;
        Tminus(i, i-1) = conj(sqrt(p*(2*j+1-p)));
        
    end

    % iterating over the elements of the (+1th) diagonal
    for i=1:(length(mvalues)-1)
        p = i;
        Tplus(i, i+1) = sqrt(p*(2*j+1-p));
    end

    i = sqrt(-1);
    
    T1 = (1/2) * (Tplus + Tminus);
    T2 = (1/(2*i)) * (Tplus - Tminus);
    T3 = diag(-mvalues);

end
