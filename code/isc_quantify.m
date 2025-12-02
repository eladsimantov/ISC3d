function [varargout] = isc_quantify(alphaThigh,alphaShank,alphaFoot,opts)
    % QuantifyISC is a function to quantify intersegmental
    % coordination. Enter the type of result you wish to quantify
    % along with the thigh shank and foot trajectories.
    
    % Usage:
    %   planarityIndex = isc_quantify(alphaThigh,alphaShank,alphaFoot,type,"PI")
    %   u3t = isc_quantify(alphaThigh,alphaShank,alphaFoot,type,"u3t")
    %   PVPC2 = isc_quantify(alphThigh,alphaShank,alphaFoot,type,"PVPC2")
    % --------------------------------------------- % 
    % It is based on definitions provided in:
    % N. A. Borghese, L. Bianchi, and Francesco Lacquaniti, 
    % “Kinematic determinants of human locomotion,” 
    % The Journal of Physiology, vol. 494, no. 3, pp. 863–879, 
    % Aug. 1996, doi: 10.1113/jphysiol.1996.sp021539.
    % AND 
    % Simon Israeli‐Korn et al., 
    % “Intersegmental coordination patterns are differently affected 
    % in Parkinson’s disease and cerebellar ataxia,” 
    % Journal of Neurophysiology, vol. 121, no. 2, pp. 672–689, 
    % Feb. 2019, doi: 10.1152/jn.00788.2017.
    % --------------------------------------------- % 
    arguments (Input)
        alphaThigh (:,1) double
        alphaShank (:,1) double
        alphaFoot (:,1) double
        opts.type {mustBeMember(opts.type, ...
            {'PI','LI','PVPC2','Eccent','U','Pearson','u3','u3t'})} = 'PI'
    end
    [V,~,latent] = pca([alphaThigh,alphaShank,alphaFoot]);
    switch opts.type
        case 'LI' 
            % linearity index
            LI = @(l1,l2,l3)(100*(l1)/(l1+l2+l3));
            varargout{1} = LI(latent(1),latent(2),latent(3));
        case 'PI' 
            % planarity index
            PI = @(l1,l2,l3)(100*(l1+l2)/(l1+l2+l3));
            varargout{1} = PI(latent(1),latent(2),latent(3));
        case 'PVPC2'
            % PVPC2 - Percentage Variance of PC2 represents the "width" of
            % the loop 
            PVPC2 = @(l1,l2,l3)(100*(l2)/(l1+l2+l3));
            varargout{1} = PVPC2(latent(1),latent(2),latent(3));
        case 'Pearson'
            % this provides the correlation between shank and foot
            % only, although you can theoretically get the
            % correlation coefficients between other angles.
            [R,P] = corrcoef(alphaFoot,alphaShank);
            varargout{1} = R(1,2);
            varargout{2} = P(1,2);       
        case 'U'
            [u1,u2,u3] = unpackMat(V);
            if u3(3) < 0, u3 = -u3; end % Correct the normal to plane to be above the plane (foot axis (=3) direction cosine is positive)
            varargout{1} = u1;
            varargout{2} = u2;
            varargout{3} = u3;
        case 'u3'
            [~,~,u3] = unpackMat(V);
            if u3(3) < 0, u3 = -u3; end % Correct the normal to plane to be above the plane (foot axis (=3) direction cosine is positive)
            varargout{1} = u3;
        case 'u3t'
            [~,~,u3] = unpackMat(V);
            if u3(3) < 0, u3 = -u3; end % Correct the normal to plane to be above the plane (foot axis (=3) direction cosine is positive)
            varargout{1} = u3(1);
        case 'Eccent'
            % Eccentricity of the ellipse - see both Barliya and Simon Israeli Korn.
            eccent = @(l1,l2)(sqrt(1-l2^2/l1^2));
            varargout{1} = eccent(latent(1),latent(2));
    end
    
    % u1t - The direction cosine of the first principle component vector
    % with respect to the thigh axis. An increase in this value
    % represents the "rotation about the u3 axis" which is observed
    % in increase in incline in Able bodied treadmill gait. 
    
    % PVPC2 - Percentage Variance of PC2 represents the "width" of
    % the loop which we have observed to be increasing with incline
    % using Able bodied treadmill data. This is in addition to the
    % rotation of the loop represented by u1t increasing which is
    % defined above.
end



function varargout = unpackMat(inputMat)
            % unpackMat Unpacks columns of a matrix into separate output variables.
            %
            %   [a, b, c, ...] = unpackMat(inputMat) assigns each column of inputMat 
            %   to a separate output variable. The number of outputs must not exceed 
            %   the number of columns in inputMat.
            %
            %   Example:
            %       A = rand(100, 3);
            %       [x, y, z] = unpackMat(A);
            for col = 1:nargout
                varargout{col} = inputMat(:,col);
            end
        end
