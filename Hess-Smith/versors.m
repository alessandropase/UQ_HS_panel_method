function [n, tau] = versors (point1, point2)

tau = point2-point1;
tau = tau/norm(tau);

n = cross([tau;0], -[0 0 1]');
n = n(1:2);