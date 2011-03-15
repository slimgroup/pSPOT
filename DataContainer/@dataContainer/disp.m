function disp(x)

% Variables pre-processing
% Distributed
if x.isdist
    dist = ['yes'...
        '\n Distr. Dimension : ' num2str(x.codist.Dimension)...
        '\n Distr. Partition : [' num2str(x.codist.Partition) ']'];        
else
    dist = 'no';
end

% Vectorized
if x.veced, vec = 'yes';
else vec = 'no'; end

% Display
fprintf('[[Data Container]] \n');
fprintf(['       Dimensions : [' num2str(size(x)) ']\n']);
fprintf(['      Permutation : [' num2str(x.perm) ']\n']);
fprintf(['       Vectorized : ' vec '\n']);
fprintf(['      Distributed : ' dist '\n']);
%  fprintf('          History :\n'); disp(x.history);