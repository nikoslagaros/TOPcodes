function [E,dE] = interpolate(x, E0, Emin, penal)
    E = Emin + x.^penal * (E0-Emin);
    dE = penal * x.^(penal-1) * (E0-Emin);
end