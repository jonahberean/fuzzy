% function scalarfieldtest(N)
% scalarfieldtest.m - Testing the definition of the scalar field in the
% matrix mechanics, to identify if it is consistent with our understanding
% of the gauge symmetry
%
% Syntax: scalarfieldtest(N)
%
% Inputs:
%    N     - coordinate matrix dimension
%           
% Outputs:
%
% Other m-files required:
% Subfunctions:
% MAT-files required: none
%
% Author: Jonah Berean-Dutcher
% email: jbd@phas.ubc.ca
% June 2021; Last revision: 27-July-2021
%------------- BEGIN CODE --------------

disp('hello world')

% % setting N
% N = 3;
    
% % define the su2 generators
% [L1, L2, L3, ~, ~] = su2generators(N);

% % Defining the three symbolic matrices
% syms Rx1 [N N];
% syms Rx2 [N N];
% syms Rx3 [N N];

% % Rx1 = rand(3);
% % Rx2 = rand(3);
% % Rx3 = rand(3);

% % Substituting in the real valued matrix elements
% x1 = (1/2*(Rx1+Rx1.')+1/(2*sqrt(-1))*(Rx1.'-Rx1));
% x2 = (1/2*(Rx2+Rx2.')+1/(2*sqrt(-1))*(Rx2.'-Rx2));
% x3 = (1/2*(Rx3+Rx3.')+1/(2*sqrt(-1))*(Rx3.'-Rx3));

% % Defining the potential for use in constructing K symbolically
% % V2 = trace(-...
% % 2j*L3*x1*x2 + 2j*L3*x2*x1 - 2*x1^2 + ... % cycle 1
% % 2*L1*x2*x2*L1 + 2*L2*x1*x1*L2 - 2*L1*x2*L1*x2 - 2*L2*x1*L2*x1 - ...
% % 2*x1*L2*L1*x2 - 2*L2*x1*x2*L1 + 2*L2*x1*L1*x2 + 2*x1*L2*x2*L1 - ...
% % 2j*L1*x2*x3 + 2j*L1*x3*x2 - 2*x2^2 + ...  % cycle 2
% % 2*L2*x3*x3*L2 + 2*L3*x2*x2*L3 - 2*L2*x3*L2*x3 - 2*L3*x2*L3*x2 - ...
% % 2*x2*L3*L2*x3 - 2*L3*x2*x3*L2 + 2*L3*x2*L2*x3 + 2*x2*L3*x3*L2 - ...
% % 2j*L2*x3*x1 + 2j*L2*x1*x3 - 2*x3^2 + ... % cycle 3
% % 2*L3*x1*x1*L3 + 2*L1*x3*x3*L1 - 2*L3*x1*L3*x1 - 2*L1*x3*L1*x3 - ...
% % 2*x3*L1*L3*x1 - 2*L1*x3*x1*L3 + 2*L1*x3*L3*x1 + 2*x3*L1*x1*L3);

% % Computing phi, the scalar field, or radial component of the matrices,
% phi = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);

% % Computing trace(phi)
% % phi = trace((1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3));

% % Computing coefficients for an expansion of phi in terms of the original
% % symbolic matrices.
% % 
% % a = (1/N)*trace(x1' * phi);
% % b = (1/N)*trace(x2' * phi);
% % c = (1/N)*trace(x3' * phi);
% % 
% % phiExpand = a*x1 + b*x2 + c*x3;

% % A = zeros(3, N, N);

% dx = zeros(3,N,N);
% fluc = Rx1;
% fluc(:,:,2) = Rx2;
% fluc(:,:,3) = Rx3;
% for i = 1:N
    
%     for j = 1:N
        
%         for k = 1:3                     
            
%             for l = 1:N
                
%                 for m = 1:N
                                        
%                     val = diff(phi(i,j), fluc(l,m,k));
%                     if 
%                     dx(k,l,m) = val;                    
                    
%                 end
%             end
%         end
%     end
% end

                
        


% % for k = 1:3
% % 
% %     for j = 1:N^2
% %         
% %         for k
% %         
% % 
% %         A(j,i) = diff(Y_phi(j), Y(i));
% % 
% %     end
% % 
% % end

% %%

% % Constructing Y, the column vector of 3N^2 length, containing the entries
% % of the three matrices.
% Y = sym(zeros(1,3*N^2));

% % Append the main diagonals to Y
% Y(1:N)       = diag(Rx1);
% Y(N+1:2*N)   = diag(Rx2);
% Y(2*N+1:3*N) = diag(Rx3);

% % index of the first element after the main diagonals
% start = 3*N+1;

% % placeholder for counting indices
% step = 0;

% % Looping from 1 to N-1, the number of diagonals of a given sign i.e. 
% % upper or lower
% for i = 1:N
%     start = start + step;
%     Y(start           : start+(N-1-i))         = diag(Rx1, i);
%     Y(start+1*(N-i)   : start+2*(N-i)-1)       = diag(Rx2, i);
%     Y(start+2*(N-i)   : start+3*(N-i)-1)       = diag(Rx3, i);
%     Y(start+3*(N-i)   : start+4*(N-i)-1)       = diag(Rx1, -i);
%     Y(start+4*(N-i)   : start+5*(N-i)-1)       = diag(Rx2, -i);
%     Y(start+5*(N-i)   : start+6*(N-i)-1)       = diag(Rx3, -i);

%     step = 6*(N-i);
% end


% % Constructing Y_phi, the column vector of N^2 length, with just the
% % entries of phi, as symbolic expressions of the entries of the three
% % matrices.
% Y_phi = sym(zeros(1,N^2));
% 
% % Append the main diagonal to Y_phi
% Y_phi(1:N)       = diag(phi);
% 
% % index of the first element after the main diagonals
% start = N+1;
% 
% % placeholder for counting indices
% step = 0;
% 
% % Looping from 1 to N-1, the number of diagonals of a given sign i.e. 
% % upper or lower
% for i = 1:N
%     start = start + step;
%     Y_phi(start           : start+(N-1-i))         = diag(phi, i);
%     Y_phi(start+1*(N-i)   : start+2*(N-i)-1)       = diag(phi, -i);
% 
%     step = 2*(N-i);
% end

% Construct the matrix that carries out the Y -> Y_phi transformation by
% taking partial derivatives
A = zeros(N^2, 3*N^2);

for i = 1:3*N^2

    for j = 1:N^2

        A(j,i) = diff(Y_phi(j), Y(i));

    end

end
%     

%     K = double((1/2)*hessian(V2, Y));

    % loading and diagonalizing K
M = load([fileparts(pwd), '/savedcouplingmatrices/',num2str(N), '.mat']);
K = real(M.K);
[V, D] = eig(K);

% retrieving the null eigenvectors
V0 = V(:, abs(diag(D)) < 1e-4);

% disp(size(V0))
%     
%     % !!! WIP, examining potential clues from symmetry in the null space of
%     % K
%     B = 1i * symmetrymatrix(N);  
%     [VB, DB] = eig(B);
%     [VK, DK] = eig(K);
%     
%     KBD = VB' * K * VB;
%     BBD = VK' * B * VK;
%     
%     [U, D] = blkbyblkdiag(BBD,DK);
% %     sym(round([diag(D), diag(DK)],1))
%     
%     
%     

% Construct T_0 from these null eigenvectors
T0 = zeros(3*N^2, 3*N^2);
for jj = 1:size(V0,2)

    T0 = T0 + V0(:,ii) * V0(:,ii)';

end

%%
r = rank(A*(eye(size(T0))-T0));
disp(r)
    
%     r = rank(A*T_0);
%     disp(r)
% 
%     
% %     
% % end

% plot(N_vec, fracNull)
% hold on
% plot(N_vec, fracRank)

% what is the size of the null space of K?
% what fraction of the total 3N^2 degrees of freedom does that represent?
% what is the rank of A*T_0?
% % how rank-deficient is that as a fraction of the potential full rank?
% 
% N = 14;
% 
% % alpha = sqrt(4 / (N*N-1));
% 
% [L1, L2, L3, ~, ~] = su2generators(N);
% 
% % L1 = alpha * L1;
% % L2 = alpha * L2;
% % L3 = alpha * L3;
% 
% % Defining the three symbolic matrices
% syms Rx1 [N N];
% syms Rx2 [N N];
% syms Rx3 [N N];
% 
% % Substituting in the real valued matrix elements
% x1 = (1/2*(Rx1+Rx1.')+1/(2*sqrt(-1))*(Rx1.'-Rx1));
% x2 = (1/2*(Rx2+Rx2.')+1/(2*sqrt(-1))*(Rx2.'-Rx2));
% x3 = (1/2*(Rx3+Rx3.')+1/(2*sqrt(-1))*(Rx3.'-Rx3));
% 
% % returns the symbolic expression for V2 in terms of the real valued 
% % entries of fluctuation matrices
% V2 = trace(-...
%     2j*L3*x1*x2 + 2j*L3*x2*x1 - 2*x1^2 + ... % cycle 1
%     2*L1*x2*x2*L1 + 2*L2*x1*x1*L2 - 2*L1*x2*L1*x2 - 2*L2*x1*L2*x1 - ...
%     2*x1*L2*L1*x2 - 2*L2*x1*x2*L1 + 2*L2*x1*L1*x2 + 2*x1*L2*x2*L1 - ...
%     2j*L1*x2*x3 + 2j*L1*x3*x2 - 2*x2^2 + ...  % cycle 2
%     2*L2*x3*x3*L2 + 2*L3*x2*x2*L3 - 2*L2*x3*L2*x3 - 2*L3*x2*L3*x2 - ...
%     2*x2*L3*L2*x3 - 2*L3*x2*x3*L2 + 2*L3*x2*L2*x3 + 2*x2*L3*x3*L2 - ...
%     2j*L2*x3*x1 + 2j*L2*x1*x3 - 2*x3^2 + ... % cycle 3
%     2*L3*x1*x1*L3 + 2*L1*x3*x3*L1 - 2*L3*x1*L3*x1 - 2*L1*x3*L1*x3 - ...
%     2*x3*L1*L3*x1 - 2*L1*x3*x1*L3 + 2*L1*x3*L3*x1 + 2*x3*L1*x1*L3);
% 
% % Computing phi, the scalar field, or radial component of the matrices,
% phi = (1/2) * (L1 * x1 + x1 * L1 + L2 * x2 + x2 * L2 + L3 * x3 + x3 * L3);
% 
% %%
% % Constructing Y, the column vector of 3N^2 length, containing the entries
% % of the three matrices.
% Y = sym(zeros(1,3*N^2));
% 
% % Append the main diagonals to Y
% Y(1:N)       = diag(Rx1);
% Y(N+1:2*N)   = diag(Rx2);
% Y(2*N+1:3*N) = diag(Rx3);
% 
% % index of the first element after the main diagonals
% start = 3*N+1;
% 
% % placeholder for counting indices
% step = 0;
% 
% % Looping from 1 to N-1, the number of diagonals of a given sign i.e. 
% % upper or lower
% for i = 1:N
%     start = start + step;
%     Y(start           : start+(N-1-i))         = diag(Rx1, i);
%     Y(start+1*(N-i)   : start+2*(N-i)-1)       = diag(Rx2, i);
%     Y(start+2*(N-i)   : start+3*(N-i)-1)       = diag(Rx3, i);
%     Y(start+3*(N-i)   : start+4*(N-i)-1)       = diag(Rx1, -i);
%     Y(start+4*(N-i)   : start+5*(N-i)-1)       = diag(Rx2, -i);
%     Y(start+5*(N-i)   : start+6*(N-i)-1)       = diag(Rx3, -i);
% 
%     step = 6*(N-i);
% end
% 
% % %%
% % % Constructing Y_const, the column vector of 3N^2 length, containing the entries
% % % of the three matrices, for the constant field zero mode
% % Y_const = zeros(1,3*N^2);
% % 
% % y = 1*randn(1,3);
% % y = bsxfun(@rdivide,y,sqrt(sum(y.^2,2)));
% % 
% % L1 = y(1) * eye(N);
% % L2 = y(2) * eye(N);
% % L3 = y(3) * eye(N);
% % 
% % L1_r = 1/2*(L1 + L1.')+1/(2*sqrt(-1))*(L1 - L1.');
% % L2_r = 1/2*(L2 + L2.')+1/(2*sqrt(-1))*(L2 - L2.');
% % L3_r = 1/2*(L3 + L3.')+1/(2*sqrt(-1))*(L3 - L3.');
% % 
% % % Append the main diagonals to Y_const
% % Y_const(1:N)       = diag(L1_r);
% % Y_const(N+1:2*N)   = diag(L2_r);
% % Y_const(2*N+1:3*N) = diag(L3_r);
% % 
% % % index of the first element after the main diagonals
% % start = 3*N+1;
% % 
% % % placeholder for counting indices
% % step = 0;
% % 
% % % Looping from 1 to N-1, the number of diagonals of a given sign i.e. 
% % % upper or lower
% % for i = 1:N
% %     start = start + step;
% %     Y_const(start           : start+(N-1-i))         = diag(L1_r, i);
% %     Y_const(start+1*(N-i)   : start+2*(N-i)-1)       = diag(L2_r, i);
% %     Y_const(start+2*(N-i)   : start+3*(N-i)-1)       = diag(L3_r, i);
% %     Y_const(start+3*(N-i)   : start+4*(N-i)-1)       = diag(L1_r, -i);
% %     Y_const(start+4*(N-i)   : start+5*(N-i)-1)       = diag(L2_r, -i);
% %     Y_const(start+5*(N-i)   : start+6*(N-i)-1)       = diag(L3_r, -i);
% % 
% %     step = 6*(N-i);
% % end
% 
% %%
% % K is the matrix for which Y.' * K * Y recovers the potential V2
% % The extra 1/2 factor is necessary. Consider the duplicate terms 
% % generated by hessian() due to symmetry of the derivatives.
% K = double((1/2)*hessian(V2, Y));
% % 


% Construct the operator T_0 which is a projection operator onto the null
% space of the coupling matrix K,
% M = load([fileparts(pwd), '/savedcouplingmatrices/',num2str(N), '.mat']);
% [V, D] = eig(real(K));
% % disp(diag(D))
% V0 = V(:, abs(diag(D)) < 1e-4);
% disp(size(V0))
% 
% % Construct T_0 from these null eigenvectors
% T_0 = zeros(3*N^2, 3*N^2);
% for i = 1:size(V0, 2)
%     T_0 = T_0 + V0(:,i) * V0(:,i)';
% end
% 
% VT = V(:, abs(diag(D)) > 1e-4);
% % Construct T_0T from these non-null eigenvectors
% T_0T = zeros(3*N^2, 3*N^2);
% for i = 1:size(VT, 2)
%     T_0T = T_0T + VT(:,i) * VT(:,i)';
% end
% 
% disp(sym(round(A*T_0, 2)))
% disp(sym(round(A, 2)))
% disp(sym(round(A*T_0T, 2)))





