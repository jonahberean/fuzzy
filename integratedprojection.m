function P = integratedprojection(N, capAngle, Z, subspaceVectors)
%{
integratedprojection.m - Returns an integrated projection operator over several projections that are interpreted to lie within a spherical cap.

Inputs:
   N          - coordinate matrix dimension
   capAngle   - angle of inclination that defines the spherical cap. It is      defined for theta in [-pi/2, pi/2] as angle made with the x-y plane. This is later changed to [pi, 0] to be in agreement with the convention for calculating spherical harmonics.

Outputs:
    P         - integrated projection operator

Other m-files required: none

Subfunctions: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
June 2021; Last revision: 05-Aug-2021 
%}

% Initializes an array to hold the integrated projection operator.
P = zeros(size(subspaceVectors,1));

% Computes the height at which to place the cut-off of the spherical cap over which we will integrate projection operators. Here we switch from the MATLAB convention to the spherical harmonic convention.
capHeight = (N-1) * cos(-capAngle + pi/2);

% Iterates over every minimal eigenvector computed, and adds the corresponding projection operator contribution to P if the height is within the cap.
for i = 1:size(Z,1)
    for j = 1:size(Z,2)
        if (Z(i,j) >= capHeight)
            P = P + subspaceVectors(:,i,j)*subspaceVectors(:,i,j)';
        end
    end
end





