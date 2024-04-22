classdef frictionTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function testZeroNetForces(testCase)
            % Test case for zero net forces
            
            netForces = [0; 0; 0]; % Zero net forces
            frictionForce = 0.5;   % Friction force magnitude
            
            expectedForces = [0; 0; 0]; % Expect zero resulting forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces);
        end
        
        function testZeroFrictionForce(testCase)
            % Test case for zero friction force
            
            netForces = [1; 2; 3]; % Example net forces
            frictionForce = 0;     % Zero friction force
            
            expectedForces = netForces; % Expect no change in forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces);
        end
        
        function testOppositeFrictionForce(testCase)
            % Test case for friction force in opposite direction
            
            netForces = [1; 0; 0]; % Example net forces
            frictionForce = -0.5;  % Friction force magnitude
            
            expectedForces = [0.5; 0; 0]; % Expect reduction in net force
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces);
        end
        
        function testFrictionGreaterThanNetForce(testCase)
            % Test case for friction force greater than net force
            
            netForces = [1; 0; 0]; % Example net forces
            frictionForce = 2;     % Friction force magnitude
            
            expectedForces = [0; 0; 0]; % Expect zero resulting forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces);
        end

                function testMultipleComponentsFriction(testCase)
            % Test case with multiple components of force and nonzero friction
            
            netForces = [3; 4; 0]; % Example net forces
            frictionForce = 2;     % Friction force magnitude
            
            expectedForces = [1; 2; 0]; % Expect reduction in net force
            
            % Call the friction function
            resultingForces = friction(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces, 'AbsTol', 1e-6);
        end
        
        function testMultipleVectorsFriction(testCase)
            % Test case with multiple force vectors and nonzero friction
            
            netForces = [3, 0, 0; 0, 4, 0; 0, 0, 5]; % Example net forces (three vectors)
            frictionForce = 2;                        % Friction force magnitude
            
            expectedForces = [1, 0, 0; 0, 2, 0; 0, 0, 3]; % Expect reduction in net forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces, 'AbsTol', 1e-6);
        end
        
        function testNegativeFrictionForce(testCase)
            % Test case with negative friction force
            
            netForces = [3; 4; 0]; % Example net forces
            frictionForce = 3;    % Negative friction force magnitude
            
            expectedForces = [1.2; 1.6; 0]; % Expect no change in forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces, 'AbsTol', 1e-6);
        end

        function testMultipleComponentVector(testCase)
            % Test case with negative friction force
            
            netForces = [2;1;0]; % Example net forces
            frictionForce = 2;    % Negative friction force magnitude
            
            expectedForces = [0.2111;0.1056;0]; % Expect no change in forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces, 'AbsTol', 1e-4);
        end

        function testMatrixEntry(testCase)
            % Test case with negative friction force
            
            netForces = [2 3;1 4;0 0]; % Example net forces
            frictionForce = 2;    % Negative friction force magnitude
            
            expectedForces = [0.2111 1.8;0.1056 2.4;0 0]; % Expect no change in forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces, 'AbsTol', 1e-4);
        end

        function testWeirdDirectionVector(testCase)
            % Test case with negative friction force
            
            netForces = [-6;3;-4]; % Example net forces
            frictionForce = 3;    % Negative friction force magnitude
            
            expectedForces = [-3.695;1.848;-2.464]; % Expect no change in forces
            
            % Call the friction function
            resultingForces = friction_3d(netForces, frictionForce);
            
            % Verify that the resulting forces match the expected forces
            verifyEqual(testCase, resultingForces, expectedForces, 'AbsTol', 1e-3);
        end

        function graphMagnitudes(testCase)
            force = [1;1;0];
            mags = -10:0.1:10;
            forces = zeros(3, length(mags));
            for i = 1:numel(mags)
                forces(:,i) = mags(i) * force;
            end
            
            figure
            plot(sqrt(sum(friction_3d(forces, 3).^2, 1)))
            hold on
            plot(sqrt(sum(forces.^2, 1)))
            hold off
            grid on
            legend('Friction magnitude', 'Frictionless Magnitude')
            title('Test to assure friction magnitude is applied correctly')
            ylabel('Magnitude')
            xlabel('Index')
        end
    end

end