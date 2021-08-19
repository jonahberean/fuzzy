function sphHarmonicMatrices = getSphericalHarmonicMatrices(N)

filename = ['savedSphericalHarmonicMatrices/',num2str(N), '.mat'];
if isfile(filename)
    M                 = load(filename);
    sphHarmonicMatrices = M.sphHarmonicMatrices;
else    

    [~, ~, ~, ~, Jplus] = computeSu2(N);

    sphHarmonicMatrices = zeros(N, N, N, (2 * N - 1));

    % Initializes an array of values for j
    jArr = 0:1:(N-1);
    
    % Iterating over j values
    for jIter = 1:length(jArr)

        % Sets the value of j
        j = jArr(jIter);

        % Computes the constant C, such that the normalization condition is satisfied.
        C = (N / (trace((Jminus^j)' * Jminus^j)))^(1/2);
        sphHarmonicMatrices(:,:,jIter,1) = C * (Jminus)^j;
        
        mArr = -j:1:j;
        
        for mIter = 1:(length(mArr)-1)

            m = mArr(mIter);
        
            % compute the Y matrices for m = -j + 1 ... j using commutators with Lplus
            sphHarmonicMatrices(:,:,jIter,mIter+1) = ...
            computeCommutator(Jplus ,sphHarmonicMatrices(:,:,jIter,mIter)) ...
            / sqrt((j - m) * (j + m + 1));

        end

    end

    % save the results
    save(filename,'sphHarmonicMatrices');

end