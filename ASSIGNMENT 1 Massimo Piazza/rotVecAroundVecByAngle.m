function v_rot = rotVecAroundVecByAngle(v,k,theta)

% PROTOTYPE:  v_rot = rotVecAroundVecByAngle(v,k,theta)
%
% DESCRIPTION:
% Rotates array of 3D vectors by an angle theta about vector k.
% Direction is determined by the right-hand rule.
%
% INPUT:
%    v - Array of three dimensional vectors to rotate. The array can be 
%           composed of N rows of 3D row vectors or N columns of 3D column
%           vectors. If v is 3x3 array, it is assumed that it is 3 rows of
%           3 3D row vectors.
%    k - Rotation axis (not necessarily a unit vector)
%    theta - Rotation angle in [rad]; positive according to right-hand
%           rule
%
%   Note: k and individual 3D vectors in v array must be same orientation.
%           
%
% OUTPUT:
%    v_rot - Array of rotated vectors.
%
%

    [m,n] = size(v);
    if (m ~= 3 && n ~= 3)
        error('input vector is/are not three dimensional')
    end
    if (size(v) ~= size(k)) 
        error('rotation vector v and axis k have different dimensions')
    end
    
    k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2); % normalize rotation axis
    No = numel(v)/3; % number of vectors in array
    v_rot = v; % initialize rotated vector array
    if ( n == 3 )
        crosskv = v(1,:); % initialize cross product k and v with right dim.
        for i = 1:No
            crosskv(1) = k(2)*v(i,3) - k(3)*v(i,2);
            crosskv(2) = k(3)*v(i,1) - k(1)*v(i,3); 
            crosskv(3) = k(1)*v(i,2) - k(2)*v(i,1);
            v_rot(i,:) = cos(theta)*v(i,:) + (crosskv)*sin(theta)...
                            + k*(dot(k,v(i,:)))*(1 - cos(theta));
        end
    else % i.e.: if m == 3 && n ~= 3
        crosskv = v(:,1); % initialize cross product k and v with right dim.
        for i = 1:No
            crosskv(1) = k(2)*v(3,i) - k(3)*v(2,i);
            crosskv(2) = k(3)*v(1,i) - k(1)*v(3,i); 
            crosskv(3) = k(1)*v(2,i) - k(2)*v(1,i);
            v_rot(:,i) = cos(theta)*v(:,i) + (crosskv)*sin(theta)...
                            + k*(dot(k,v(:,i)))*(1 - cos(theta));
        end
    end
end