function phi = calculate_phase(T, t_ph, t_pQ)
    % CALCULATE_PHASE Computes the phase difference
    % 
    % Inputs:
    %   T    - Period of the oscillation
    %   t_ph - Time corresponding to phase measurement
    %   t_pQ - Reference time for phase calculation
    %
    % Output:
    %   phi  - Calculated phase difference in radians
    %
    % Formula: phi = (2 * pi / T) * (t_ph - t_pQ)

    % Compute phase difference
    phi = (2 * pi / T) * (t_ph - t_pQ);
end
