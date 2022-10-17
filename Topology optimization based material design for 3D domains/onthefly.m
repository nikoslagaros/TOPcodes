classdef onthefly < handle
    properties
        Nb % number of reduced basis vectors
        fi % reduced basis matrix
        tol % accuracy tolerance
        loop % iteration
        fll % number of full iterations
        rdc % number of reduced iterations
    end

    methods
        function obj = onthefly(nb,tol)
            obj.Nb = nb;
            obj.loop = 0;
            obj.fll = 0;
            obj.rdc = 0;
            obj.tol = tol;
        end

        function U = solve(obj,K,F)
            obj.loop = obj.loop + 1;
            if obj.loop == 1
                U = obj.fea(K,F);
                obj.fi(:,obj.loop) = U/norm(U);
            elseif obj.loop <= obj.Nb
                U = obj.fea(K,F);
                Uorth = U - obj.fi*(obj.fi'*U);
                obj.fi(:,obj.loop) = Uorth/norm(Uorth);
            else
                y = obj.fi'*K*obj.fi \ obj.fi'*F;
                U = obj.fi*y;
                dF = K*U-F;
                res = norm(dF);
                if res > obj.tol
                    U = obj.fea(K,F);
                    obj.fi(:,1) = [];
                    Uorth = U - obj.fi*(obj.fi'*U);
                    obj.fi(:,obj.Nb) = Uorth/norm(Uorth);
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

                function flag = sameSize(obj, sz)
            flag = size(obj.fi,1) == sz;
        end

        function p = reinitialize(obj)
            p = onthefly(obj.Nb, obj.tol);
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