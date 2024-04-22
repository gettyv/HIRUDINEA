function resultingForces = friction(netForces, frictionForces)
%FRICTION Apply a friction force or torque to the given net force or torque
%   Applies the friction in the opposite direction and ensures
%   that the friction magnitude does not exceed the magnitude of the
%   combination for each element in the arrays

% Calculate the minimum magnitude between the absolute values of netForces
% and frictionForces
clampedFrictionMag = min(abs(netForces), abs(frictionForces));

% Apply friction based on the sign of netForces
combined = netForces - sign(netForces) .* clampedFrictionMag;

% Store the resulting forces
resultingForces = combined;

end