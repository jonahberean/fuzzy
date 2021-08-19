function field = getFieldFromCoefficients(c, sphHarmonics)

field = zeros(size(sphHarmonics, 1), size(sphHarmonics, 2));

for j = 0:(length(c)-1)    
    
    mValues = -j:1:j;
    
    for m = 1:length(mValues)
        
        for p = 1:size(sphHarmonics, 2)

            field = field + c(j+1,m) * sphHarmonics(:, :, j, m)
            
        end
        
    end
    
end