classdef pod < handle
    properties
        Nb % number of reduced basis vectors
        fi % reduced basis matrix
        A  % displacement snapshots
        tol % accuracy tolerance
        loop % iteration
        fll % number of full iterations
        rdc % number of reduced iterations
    end

    methods
        function obj = pod(nb,tol)
            obj.Nb = nb;
            obj.loop = 0;
            obj.fll = 0;
            obj.rdc = 0;
            obj.tol = tol;
        end

        function U = solve(obj,K,F)
            obj.loop = obj.loop + 1;
            if obj.loop <= obj.Nb
                U = obj.fea(K,F);
                obj.A(:,obj.loop) = U;
                if obj.loop == obj.Nb
                    obj.A = obj.A - mean(U);
                    [obj.fi,~,~] = svd(obj.A,'econ');
                end
            else
                y = obj.fi'*K*obj.fi \ obj.fi'*F;
                U = obj.fi*y;
                dF = K*U-F;
                res = norm(dF);
                if res > obj.tol
                    U = obj.fea(K,F);
                    obj.A(:,1) = [];
                    obj.A(:,obj.Nb) = U;
                    obj.A = obj.A - mean(U);
                    [obj.fi,~,~] = svd(obj.A,'econ');
                else
                    obj.rdc = obj.rdc + 1;
                end
            end

        end

        function U = fea(obj,K,F)
                fprintf('full iteration \n') 
                dK = decomposition(K);
                U = dK\F;
                obj.fll = obj.fll + 1;
        end

        function [fll,rdc] = counts(obj)
            fll = obj.fll;
            rdc = obj.rdc;
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by G. Kazakis and N.D. Lagaros   
%
% The code is intended for educational purposes, extensions can be found in
% the paper "A simple Matlab code for material design optimization using
% reduced order models"
%
% Disclaimer:                                                              
% The authors reserves all rights for the program.    
% The code may be distributed and used for educational purposes.       
% The authors do not guarantee that the code is free from errors, and
% they shall not be liable in any event caused by the use of the program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%