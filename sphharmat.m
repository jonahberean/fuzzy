function [Y, Jplus, Jminus, Ji] = sphharmat(N)
% sphharmat.m - generates the 2l+1 matrices which
% are counterparts to the spherical harmonic functions
%
% Syntax: Y = sphharmat(N)
%
% Inputs:
%    N - matrix dimension, must be odd (WHY?)
%
% Outputs:
%    Y      - spherical harmonic matrices
%    Lplus  - The raising operator for this dim. irrep
%    Lminus - The lowering operator for this dim. irrep
%
% Other m-files required:
% Subfunctions:
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% June 2021; Last revision: 01-June-2021
%------------- BEGIN CODE --------------

% Initialize arrays that hold the various testing values generated 
% JminusTest = 0;
% JplusTest = 0;
% J3Test = 0;
% laplacianTest = 0;
 
% Initializes a multi-dimensional array to hold the spherical harmonic matrix counterparts.
% 1. first matrix dimension
% 2. second matrix dimension
% 3. j index
% 4. m index
Y = zeros(N, N, N, (2 * N - 1));
 
% Returns the generators and ladder operators of the su(2) Lie algebra
[J1,J2,J3,Jminus,Jplus] = su2generators(N);

% % Displays them for testing purposes. 
% Jminus
% Jplus

% Initializes and fills a multi-dimensional array with the three generators.
Ji = zeros(size(J1,1), size(J1,2), 3);
Ji(:,:,1) = J1;
Ji(:,:,2) = J2;
Ji(:,:,3) = J3;

% Initializes an array of values for j
jArr = 0:1:(N-1);
 
% Iterating over j values
for jIter = 1:length(jArr)

    % Sets the value of j
    j = jArr(jIter);
    % fprintf('j = %d\n', j)

    % Computes the constant C, such that the normalization condition is satisfied.
    C = (N / (trace((Jminus^j)' * Jminus^j)))^(1/2);
    Y(:,:,jIter,1) = C * (Jminus)^j;

    % Testing that the normalization condition is satisfied.
    % fprintf('Normalized: (1/N)*trace(Y_{j,-j} * Y{j,-j}) = %.2f\nC = %.2f\n', (1/N)*trace(Y(:,:,jIter,1)' * Y(:,:,jIter,1)), C)
    
    % Testing [Jminus, Y_{j,(-j)}] = 0
    % JminusTest = max(abs(commutator(Jminus,Y(:,:,jIter,1))), [], 'all');
    % if (JminusTest > 0.001)
    %     statement = ...
    %     '\nWarning: Largest value in [Jminus, Y_{j,(-j)}] is %.6f \n';
    %     fprintf(statement, max(JminusTest,[],'all'))
    % else
    %     fprintf('Jminus Test has passed.\n')
    % end
    
    
    mArr = -j:1:j;
    
    for mIter = 1:(length(mArr)-1)

        m = mArr(mIter);

        % fprintf("m = %d\n", mArr(mIter))
        % fprintf("Now computing Y_{%d,%d} from Y_{%d,%d}\n", j, m+1, j, m)
    
        % compute the Y matrices for m = -j + 1 ... j using commutators with Lplus
        Y(:,:,jIter,mIter+1) = ...
        commutator(Jplus ,Y(:,:,jIter,mIter)) ...
         / sqrt((j - m) * (j + m + 1));
        




%         % Testing \sum_i [Ji, [Ji, Y_{j,m}] = j(j+1) Y_{j,m}
%         laplacianTest = zeros(size(J1));
%         for i = 1:3
            
%             com = (Ji(:,:,i) * Y(:,:,j+1,m+1) - Y(:,:,j+1,m+1) * Ji(:,:,i));
%             laplacianTest = laplacianTest + (Ji(:,:,i) * com - com * Ji(:,:,i));
            
%         end        
%         [ii,~,v] = find(laplacianTest./ Y(:,:,j+1,m+1));
%         ii = (ii(~isnan(v)));
%         v  = (v(~isnan(v)));
%         ii = (ii(~isinf(v)));        
%         v  = (v(~isinf(v)));
%         out = accumarray(ii,v,[],@mean);
%         if (abs(out - (j*(j+1))) > 0.001)
%             statement = '\nWarning: Laplacian eigenvalue deviates from j(j+1) by %.6f\n';
% %             fprintf(statement, abs(out - (j*(j+1)))) 
%         end
                
%         % Testing [J3, Y_{j,m}] = m Y_{j,m}
%         J3Test = J3 * Y(:,:,j+1,m+1) - Y(:,:,j+1,m+1) * J3;
%         [ii,~,v] = find(J3Test./ Y(:,:,j+1,m+1));
%         ii = (ii(~isnan(v)));
%         v  = (v(~isnan(v)));
%         ii = (ii(~isinf(v)));        
%         v  = (v(~isinf(v)));
%         out = accumarray(ii,v,[],@mean);
%         if (abs(out - (mArr(m+1))) > 0.001)
%             statement = '\nWarning: J3Test deviates from m by %.6f\n';
% %             fprintf(statement, abs(out - (mArr(m+1)))) 
%         end
        
        
%         % Testing [Jplus, Y_{j,j}] = 0
%         if (m == (length(mArr)-1))
            
%             JplusTest = max(abs((Jplus * Y(:,:,j+1,m+1) - Y(:,:,j+1,m+1) * Jplus)), [], 'all');
%             if (JplusTest > 0.001)
%                 statement = '\nWarning: Largest value in [Jplus, Y_{j,j}] is %.6f \n';
% %                 fprintf(statement, max(JplusTest,[],'all'))
%             end            
            
%         end        
        
    end
    
end
 
% % Testing orthonormality of the Y_{j,m}, 
% % i.e. (1/N)*tr{Y_{j,m}' * Y_{k,n})} = delta_{j,k} * delta_{m, n}

% % !!! Orthonormality is also failing, we need to look into how bad it is
% for j = 0:jmax
    
%     mArr = -j:1:j;
    
%     for m = 1:length(mArr)
        
%         for k = 0:jmax
            
%             for n = 1:length(mArr)
                
%                 delta = (1/N)*trace(Y(:,:,j+1,m)'*Y(:,:,k+1,n));

%                 if (j == k)
                    
%                     if (m == n)
                        
%                         if (abs(delta - 1) > 0.001)
                            
%                             statement = '\nWarning: Orthonormality has failed for in {%d, %d}, {%d, %d} \n';
%                             fprintf(statement, j, m, k, n)                                                                                   
                            
%                         end                                                
                        
%                     else
                    
%                         if (abs(delta) > 0.001) 

%                             statement = '\nWarning: Orthonormality has failed for in {%d, %d}, {%d, %d} \n';
%                             fprintf(statement, j, m, k, n)    

%                         end
                    
%                     end
                    
%                 else
                    
%                     if (abs(delta) > 0.001) 

%                         statement = '\nWarning: Orthonormality has failed for in {%d, %d}, {%d, %d} \n';
%                         fprintf(statement, j, m, k, n)    

%                     end                    
                    
%                 end
                
                              
%             end
            
%         end
        
%     end
    
% end      
    
    
% end

filename = ['savedsphericalharmonicmatrices/',num2str(N), '.mat'];
% save the results
save(filename,'Y');
