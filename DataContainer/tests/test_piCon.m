function test_suite = test_piCon
initTestSuite;
end

function test_piCon_reshape
%% Testing piCon reshape
n1 = randi(5);
n2 = randi(5);
n3 = randi(5);
x = piCon.randn(n1,n2,n3);
x = reshape(x,n1,n2*n3);
x = reform(x);
end % reshape

function test_piCon_redistribute
%% Testing piCon reshape
n1 = randi(5);
n2 = randi(5);
n3 = randi(5);
x = piCon.randn(n1,n2,n3);
x.codistInfo;
x = redistribute(x,2);
x.codistInfo;
end % reshape