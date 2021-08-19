% function harmonicYtest(N)

%region - doc
%{
harmonicYtest.m - Some tests to investigate harmonicY.m's behaviour

Inputs:
   N         - size of the coordinate matrices

Outputs:

Other m-files required: harmonicY.m

Subfunctions: none
MAT-files required: none

Author: Jonah Berean-Dutcher
email: jbd@phas.ubc.ca
Aug 2021; Last revision: 18-Aug-2021 
%}
%endregion - doc

N = 8;
numFieldGrid = 15;
[x,y,z]= sphere(numFieldGrid-1);
[phi,theta,r] = cart2sph(x,y,z);
theta = theta + pi/2;

f = zeros(length(theta), length(phi));

A = complex(rand(N, N), rand(N,N));
A = (1/2) * (A + A');
[field, theta, phi, c] = matrixtofield(A, numFieldGrid);

% Computing the functioin values
for j = 0:(N-1)    

    % fprintf('Iterating over j. Currently at j = %.2f (out of %.2f).\n', j, jmax);
    
    mValues = -j:1:j;
    
    for m = 1:length(mValues)
        
        
        % Note: The function harmonicY does not evaluate a grid all at once, it evaluates f(theta, phi) once for each value of theta and of phi, as ordered pairs running down each column in tandem.

        for p = 1:length(phi)

            % Adds a column of values into f, that's at every theta value, and at one phi value. 
            % f(:,p) = f(:,p) + c(j+1,m) * ...
            %     harmonicY(j, mValues(m), ...
            %     theta, ones(1,length(theta))*phi(p), 'type', 'real')';

            f(:,p) = f(:,p) +  harmonicY(j, mValues(m), theta(:,p), phi(:,p));
        end
        
    end
    
end

c1 = c(2,1);
c2 = c(2,3);

Y1 = harmonicY(1,-1,(pi/4)*ones(1,length(phi(3,:))), phi(3,:));
Y2 = harmonicY(1,1,(pi/4)*ones(1,length(phi(3,:))), phi(3,:));
