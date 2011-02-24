% testing window functions
clear;
M=[10 11];
P=2;
H=1;

disp('start test');

for m=M(1):M(2)
    for p=2:P
        for h=0:H
    
	    try
                [A a1 a2]=pSPOT.pWindow.fd_window_1Dfor(m,p,h);
                [B b1 b2]=pSPOT.pWindow.fd_window_1Dbck(m,p,h);
                try
                    n=m+(p-1)*2*h;
                    x0=rand(m,1);
                    y=A*x0;
                    x1=B*y;
                    check=sum(abs((x1-x0)'));
                    fprintf('\tFD: m=%d n=%d p=%d h=%d: %d\n',m,n,p,h,check);
                    assert(check<1e-10,'ERROR: does not pass inverse test check=%f',check);
                catch xy
                    fprintf('FD: m=%d n=%d p=%d h=%d\n',m,n,p,h);
                    fprintf('FD: a1=%d a2=%d b1=%d b2=%d\n',sum(a1),sum(a2),sum(b1),sum(b2));
                    disp(xy.message);
                    whos a1 a2 A b1 b2 B
                    disp(full(A))
                    disp(full(B))
                    disp(full(B*A))
                end
	    catch AB
            fprintf('FD: m=%d n=%d p=%d h=%d\n',m,n,p,h);
	    	disp(AB.message);
	    end
    
	    try
                [A a1 a2]=pSPOT.pWindow.tpr_window_1Dfor(m,p,h);
                [B b1 b2]=pSPOT.pWindow.tpr_window_1Dbck(m,p,h);
                try
                    n=m+(p-1)*2*h;
                    x0=rand(m,1);
                    y=A*x0;
                    x1=B*y;
                    check=sum(abs((x1-x0)'));
                    fprintf('\tTPR: m=%d n=%d p=%d h=%d: %d\n',m,n,p,h,check);
                    assert(check<1e-10,'ERROR: does not pass inverse test check=%f',check);
                catch xy
                    fprintf('TPR: m=%d n=%d p=%d h=%d\n',m,n,p,h);
                    fprintf('TPR: a1=%d a2=%d b1=%d b2=%d\n',sum(a1),sum(a2),sum(b1),sum(b2));
                    disp(xy.message);
                    whos a1 a2 A b1 b2 B
                    disp(full(A))
                    disp(full(B))
                    disp(full(B*A))
              end
	    catch AB
            fprintf('TPR: m=%d n=%d p=%d h=%d\n',m,n,p,h);
	    	disp(AB.message);
	    end

        end
    end
end

disp('end test');
