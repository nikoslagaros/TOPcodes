function C = elactMat(E,nu)
%elemMat compute elasticity matrix
%   E: Young Modulus
%   nu: Poisson's ratio
A = [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
C = E/(1-nu^2)*A;
end