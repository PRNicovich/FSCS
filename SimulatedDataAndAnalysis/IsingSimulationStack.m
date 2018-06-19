
% Generate binary map based on square lattice Ising model
% Calculations based on ising.m from https://www.physics.ohio-state.edu/~braaten/statphys/Ising_MatLab.pdf

function IsingImgStack = IsingSimulationStack(IsingParams)



J = IsingParams.InteractionStrength; 
kB = 1.3806e-23; % Boltzmann constant
T = IsingParams.Temperature;

IsingImgStack = zeros(IsingParams.RegionSize(1), IsingParams.RegionSize(2), IsingParams.Nimages);

%% Generate a random initial configuration
grid = (rand(IsingParams.RegionSize(1), IsingParams.RegionSize(2)) > 0.5)*2 - 1;
flip = zeros(IsingParams.RegionSize(1), IsingParams.RegionSize(2)); 
flip(1:2:IsingParams.RegionSize(1), 1:2:IsingParams.RegionSize(2)) = 1.0; 
flip(2:2:IsingParams.RegionSize(1), 2:2:IsingParams.RegionSize(2)) = 1.0;

for k = 1:IsingParams.PreRunIterations
    [grid, flip] = evolveGrid(grid, flip, J, T, kB, IsingParams.RegionSize(1), IsingParams.RegionSize(2));
end

for k = 1:IsingParams.Nimages
    [grid, flip] = evolveGrid(grid, flip, J, T, kB, IsingParams.RegionSize(1), IsingParams.RegionSize(2));
    IsingImgStack(:,:,k) = (grid + 1)/2;
end





end

    function [grid, flip] = evolveGrid(grid, flip, J, T, kB, r1, r2)

        % Calculate the number of neighbors of each cell
        neighbors = circshift(grid, [ 0 1]) + ... 
            circshift(grid, [ 0 -1]) + ... 
            circshift(grid, [ 1 0]) + ... 
            circshift(grid, [-1 0]);
        % Calculate the change in energy of flipping a spin
        DeltaE = 2 * J * (grid .* neighbors);
        % Calculate the transition probabilities
        p_trans = exp(-DeltaE/(kB * T));
        % Decide which transitions will occur
        transitions = (rand(r1, r2) < p_trans ).*flip*-2 + 1;
        flip = ~flip;
        % Perform the transitions
        grid = grid .* transitions;

    end