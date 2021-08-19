N = 3;
xi = rand(9);
xi2 = rand(9);
gridSize = 5;
[varZ, Z, localStateVectors] = getLocalStatesConOpt(N, xi, xi2, gridSize);