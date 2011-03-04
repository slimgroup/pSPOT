function y = gatherchk(tmpy,mode,gather)
% GATHERCHK     Function for checking and processing of the gathering 
%               of the result
%
%   GATHER specifies whether to gather the results to a local array
%   or leave them distributed, default is 0.
%   GATHER = 0 will leave them distributed.
%   GATHER = 1 will gather the results of forwards or adjoint multiplication.
%   GATHER = 2 will gather only in forward mode.
%   GATHER = 3 will gather only in backward (adjoint) mode.

if mode == 1
    if op.gather == 1 || op.gather == 2
        y = gather(tmpy);
    else
        y = tmpy;
    end
else % mode == 2
    if op.gather == 1 || op.gather == 3
        y = gather(tmpy);
    else
        y = tmpy;
    end
end % gather