function [E,dE] = interpolate(x, E0, Emin, penal)
%INTERPOLATE MATERIAL PROPERTIES USING THE SIMP APPROACH
E = Emin + x.^penal * (E0-Emin);
dE = penal * x.^(penal-1) * (E0-Emin);
end