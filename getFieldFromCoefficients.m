function field = getFieldFromCoefficients(c, sphHarmonics)

field = zeros(size(sphHarmonics, 1), size(sphHarmonics, 2));

for j = 0:(size(c,1)-1)    
    
    mValues = -j:1:j;
    
    for m = 1:length(mValues)
        
        field = field + c(j+1,m) * sphHarmonics(:, :, j+1, m);       
        
    end
    
end