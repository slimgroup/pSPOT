% testing window functions
clear all;
%N=[30 1+randi(10) 500]; % range of vector sizes
N=[20 1 100]; % range of vector sizes
P=10;         % range of processors
H=5;         % range of half-overlap
T=13;        % tolerrance 1e-T

disp('start test');
fprintf('n=[%d:%d:%d] p=[2:%d] h=[0:%d] T=1e-%d\n',N(1),N(2),N(3),P,H,T);

for n=N(1):N(2):N(3)
    for p=2:P
        for h=0:H
    
            try
                [ m ysf xsf ] = pSPOT.pWindow.funWindowShape1D( n, p, h );
                %fprintf('\ttest: n=%d p=%d h=%d\n',n,p,h);
            catch pr
                %fprintf('shapes skipped: n=%d p=%d h=%d\n',n,p,h);
                %disp(pr.message);
                continue;
            end

            try
                A=oplWindow1Dfd(n,p,h);
                B=A';
                try
                    x0=rand(n,1);
                    y=A*x0;
                    x1=B*y;
                    check=sum(abs((x1-x0)'));
                    assert(check<10^-T,'ERROR: does not pass inverse test check=%f',check);
                catch xy
                    fprintf('FD: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                    disp(xy.message);
                    %disp(full(A))
                    %disp(full(B))
                    %disp(full(B*A))
                end
            catch AB
                fprintf('FD: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                disp(AB.message);
            end
    
            try
                A=oplWindow1Dtpr(m,p,h);
                B=A';
                try
                    x0=rand(n,1);
                    y=A*x0;
                    x1=B*y;
                    check=sum(abs((x1-x0)'));
                    assert(check<10^-T,'ERROR: does not pass inverse test check=%f',check);
                catch xy
                    fprintf('TPR: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                    disp(xy.message);
                    %disp(full(A))
                    %disp(full(B))
                    %disp(full(B*A))
                end
            catch AB
                fprintf('TPR: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                disp(AB.message);
            end

        end
    end
end

disp('end test');
