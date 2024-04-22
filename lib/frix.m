function resultingForces = frix(netForce, frictionForce)
%FRICTION Apply a friction force or torque to the given net force or torque
%   Applies the friction in the opposite direction and ensures
%   that the friction magnitude does not exceed the magnitude of the
%   combination for each vector in the arrays
if ~isequal([3,1],size(netForce))
    error('netForce vector is of incorrect size');
end
% Calculate the magnitude of the net forces
netForceMagnitude = sqrt(sum(netForce.^2));

% Calculate the minimum magnitude between the absolute values of net forces
% and the scalar friction force
clampedFrictionMag = min(netForceMagnitude, frictionForce);

% Calculate the unit vectors of the net forces
unitNetForce = netForce ./ netForceMagnitude;

% Calculate the combined forces
combined = netForce - unitNetForce .* clampedFrictionMag;

% Replace any NaN values with 0
combined(isnan(combined)) = 0;

% Store the resulting forces
resultingForces = combined;

end
