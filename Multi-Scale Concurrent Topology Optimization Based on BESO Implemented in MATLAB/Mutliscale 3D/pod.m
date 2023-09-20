classdef pod < handle
    properties
        Nb % number of reduced basis vectors
        fi % reduced basis matrix
        A  % displacement snapshots
        tol % accuracy tolerance
        loop % iteration
        full % number of full iterations
        reduced % number of reduced iterations
    end

    methods
        function obj = pod(nb,tol)
            obj.Nb = nb;
            obj.loop = 0;
            obj.full = 0;
            obj.reduced = 0;
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
                    obj.reduced = obj.reduced + 1;
                end
            end

        end

        function U = fea(obj,K,F)
                fprintf('full iteration \n') 
                dK = decomposition(K);
                U = dK\F;
                obj.full = obj.full + 1;
        end

        function [full,reduced] = iterations(obj)
            full = obj.full;
            reduced = obj.reduced;
        end

    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was modified by G. Kazakis and N.D. Lagaros and was 
% based on the pod.m code written by G. Kazakis and N.D. Lagaros, 
% in the paper "A simple Matlab code for material design optimization using
% reduced order models", G. Kazakis and N.D. Lagaros, 
% Materials, 2022   
%
% The code is intended for educational purposes, extensions can be found in
% the paper "Topology optimization based material design for 3D domains
% using MATLAB"
%
% Disclaimer:                                                              
% The authors reserves all rights for the program.    
% The code may be distributed and used for educational purposes.       
% The authors do not guarantee that the code is free from errors, and
% they shall not be liable in any event caused by the use of the program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%