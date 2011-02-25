function test_suite = test_oppFunComposite
%test_oppCompositeFun  Unit tests for the oppCompositeFun operator
initTestSuite;
end

function test_oppFunComposite_cell
%% Testing composite S with cells
spmd
    for k=1:labindex
        S{k} = randn(2*k,k);
    end
end

end %