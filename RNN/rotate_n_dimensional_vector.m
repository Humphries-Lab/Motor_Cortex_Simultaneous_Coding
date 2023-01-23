function v_prime = rotate_n_dimensional_vector(v,a)

% ROTATE_N_DIMENSIONAL_VECTOR rotates an n-D vector by a given angle
% VPRIME = ROTATE_N_DIMENSIONAL_VECTOR(V,A) rotates the n-dimensional
% vector in V by the angle A (in radians)
% Returns: VPRIME, the rotated vector
% 
% Notes:
% this approach rotates the vector in each 2-D hyperplane independently, 
% by considering rotating between the unit vectors that define each
% hyperplane. Which means:
%   (a) the rotations happens between pairwise adjacent entries in V:
%   (1,2), then (3,4) etc
%   (b) the rotation is only complete if n is even
%
% References:
%   https://analyticphysics.com/Higher%20Dimensions/Rotations%20in%20Higher%20Dimensions.htm
%
% Initial version: 21/7/2021
% Mark Humphries

% make it a column vector
[~,c] = size(v);
if c > 1
    v = v';
end

% check dimensions of vector
n = length(v);
if rem(n,2)
    warning('Rotations will not include final dimension because N is odd')
end

% rotate
v_prime = v;
pairs = reshape([1:n],2,n/2)';  % create a list of all pairs
for i = 1:floor(n/2)        % skip last dimension if n is odd
    % create unit vectors defining hyperplane between adjacent dimensions
    unit_vector1 = zeros(n,1);
    unit_vector1(pairs(i,1)) = 1;
    unit_vector2 = zeros(n,1);
    unit_vector2(pairs(i,2)) = 1;

    % sanity check they are orthogonal
%     unit_vector1'*unit_vector1 == 1;
%     unit_vector2'*unit_vector2 == 1;
%     unit_vector1'*unit_vector2 == 0;
%     unit_vector2'*unit_vector1 == 0;

    % create rotation matrix
    outerproduct1_1 = unit_vector1 * unit_vector1';
    outerproduct1_2 = unit_vector1 * unit_vector2';
    outerproduct2_1 = unit_vector2 * unit_vector1';
    outerproduct2_2 = unit_vector2 * unit_vector2';
    
    % this general rotation equation from the above reference
    % defines 2-D rotation matrix on given hyperplane
    n_d_rotation_matrix = eye(n) + (outerproduct2_1 - outerproduct1_2)*sin(a) + ...
                (outerproduct1_1 + outerproduct2_2)*(cos(a)-1);

    % rotate vector on this hyperplane
    v_prime = n_d_rotation_matrix * v_prime;
  
end

% for internal checking
cosine_similarity = v'*v_prime / (norm(v) * norm(v_prime));

