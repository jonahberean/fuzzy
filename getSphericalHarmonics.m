function sphHarmonics = getSphericalHarmonics(N, polarGrid, azimuthGrid, overwriteFlag)

filename = ['savedSphericalHarmonics/',num2str(N), '.mat'];
if isfile(filename) & ~strcmp('overwrite', overwriteFlag)
    M                 = load(filename);
    sphHarmonics = M.sphHarmonics;
else        
    sphHarmonics = zeros(length(polarGrid), length(azimuthGrid), N, (2 * N - 1));

    for j = 0:(N-1)    

        fprintf('Generating spherical harmonics; j = %d\n', j)
        
        mValues = -j:1:j;
        
        for m = 1:length(mValues)
            
            for p = 1:length(azimuthGrid)

                sphHarmonics(: , p, j+1, m) =  ...
                computeSphericalHarmonic(j, mValues(m), ...
                polarGrid, azimuthGrid(p) * ones(size(polarGrid)));
                
            end
            
        end
        
    end

    % save the results
    save(filename,'sphHarmonics');
    
end
