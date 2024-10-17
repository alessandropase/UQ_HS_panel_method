function output = ROT_alpha(array, alpha)

% INPUT: - array: matrix Nx2. Every row contain the coordinates of one
%                 point of the foil
%        - alpha: AoA of the foil (positive nose up)    [rad]
% 
% OUTPUT: matrix Nx2. Every row contain the coordinates of one
%         point of the foil, rotated by alpha


ROT_matrix = [cos(alpha), sin(alpha);
            -sin(alpha), cos(alpha)];
xx = array(:, 1);
yy = array(:, 2);

N = length(xx);

output = zeros(size(array));

for i = 1 : N
    vect = [xx(i); yy(i)];
    vect_rot = ROT_matrix * vect;
    output(i, 1) = vect_rot(1);
    output(i, 2) = vect_rot(2);
end
