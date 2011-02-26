function test_suite = test_oppFunComposite
%test_oppCompositeFun  Unit tests for the oppCompositeFun operator
initTestSuite;
end

function test_oppFunComposite_cell
%% Testing composite S with cells
% Multiple vectors of different dimensions stored in the cells stored as a 
% composite. Wonderful

spmd % spmd way of generating composites
    S1 = cell(1,labindex);
    for k=1:labindex
        S1{k} = randn(2*k,k);
    end
end
disp(S1);

% Local way of generating composites
% No spmd!!! HUZZAH!!!
A1 = {randn(2,1)};
A2 = {randn(4,1) randn(4,1)};
S2 = composite;
S2{1} = A1;
S2{2} = A2;
disp(S2);
end %