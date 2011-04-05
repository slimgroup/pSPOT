function test_suite = test_piCon
initTestSuite;
end

function test_piCon_reshape
%% Testing piCon reshape
x = piCon.randn(randi(10),randi(10),randi(10));
x = vec(x);
x = reform(x);
end % reshape