function resultingForces = friction_3d(netForces, frictionForce)
%FRICTION Apply a friction force or torque to the given net force or torque
%   Applies the friction in the opposite direction and ensures
%   that the friction magnitude does not exceed the magnitude of the
%   combination for each vector in the arrays

% Calculate the magnitude of the net forces
netForceMagnitudes = sqrt(sum(netForces.^2, 1));

% Calculate the minimum magnitude between the absolute values of net forces
% and the scalar friction force
clampedFrictionMag = min(netForceMagnitudes, abs(frictionForce));

% Calculate the unit vectors of the net forces
unitNetForces = netForces ./ netForceMagnitudes;

% Calculate the combined forces
combined = netForces - unitNetForces .* clampedFrictionMag;

% Replace any NaN values with 0
combined(isnan(combined)) = 0;

% Store the resulting forces
resultingForces = combined;

end
