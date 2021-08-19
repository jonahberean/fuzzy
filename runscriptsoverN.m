function runscriptsoverN(startN, stopN)

%region - doc
%{
runscriptsoverN.m - 

Inputs:


Outputs:


Other m-files required: none
Subfunctions: none
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
Aug 2021; Last revision: 18-Aug-2021 
%}
%endregion - doc

% Defines an array of N values and iterates over them
NArr = startN:1:stopN;
for i = 1:length(NArr)

    N = NArr(i);
    fprintf('Now at N = %d... going to %d.\n', N, stopN)
    

    % Generate coupling matrices symbolically.
    % K = couplingmatrixsymbolic(N);

    % Generate scalar field projection matrices
    % A = scalarfieldprojection(N);

    % Generate spherical harmonic matrices
    fprintf('Generating spherical harmonic matrices.\n');
    [Y, Jplus, Jminus, Ji] = sphharmat(N);


end