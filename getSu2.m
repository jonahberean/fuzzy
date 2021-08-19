function [J1, J2, J3,Jminus,Jplus] = getSu2(N)

j = (N-1)/2;

mvalues = -j:1:j;

Jminus = zeros(length(mvalues));
Jplus  = zeros(length(mvalues));

for i=2:(length(mvalues))
    p = i-1;
    Jminus(i, i-1) = conj(sqrt(p*(2*j+1-p)));
    
end

for i=1:(length(mvalues)-1)
    p = i;
    Jplus(i, i+1) = sqrt(p*(2*j+1-p));
end

i = sqrt(-1);

J1 = (1/2) * (Jplus + Jminus);
J2 = (1/(2*i)) * (Jplus - Jminus);
J3 = diag(-mvalues);