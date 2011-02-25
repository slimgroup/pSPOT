% testing window functions
clear all;
M=[30 10 150]; % range of vector sizes
P=4;            % range of processors
H=3;            % range of half-overlap
T=13;          % tolerrance 1e-T

disp('start test');
fprintf('n=[%d:%d:%d] p=[2:%d] h=[0:%d] T=1e-%d\n',M(1),M(2),M(3),P,H,T);

for m=M(1):M(2):M(3)
    for p=2:P
        for h=0:H
    
	    try
	        [ n w0 w1 w2 r0 r d0 ] = pSPOT.pWindow.window1Dparams( m, p, h );
	    catch pr
		fprintf('FD: m=%d p=%d h=%d\n',m,p,h);
	        disp(pr.message);
		continue;
	    end

	    try
                [A a1 a2]=pSPOT.pWindow.fd_window_1Dfor(m,p,h);
                [B b1 b2]=pSPOT.pWindow.fd_window_1Dbck(m,p,h);
                try
                    n=m+(p-1)*2*h;
                    x0=rand(m,1);
                    y=A*x0;
                    x1=B*y;
                    check=sum(abs((x1-x0)'));
                    %fprintf('\tFD: m=%d n=%d p=%d h=%d: %d\n',m,n,p,h,check);
                    assert(check<10^-T,'ERROR: does not pass inverse test check=%f',check);
                catch xy
                    fprintf('FD: m=%d n=%d p=%d h=%d\n',m,n,p,h);
                    fprintf('FD: a1=%d a2=%d b1=%d b2=%d\n',sum(a1),sum(a2),sum(b1),sum(b2));
                    disp(xy.message);
                    whos a1 a2 A b1 b2 B
                    %disp(full(A))
                    %disp(full(B))
                    %disp(full(B*A))
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
                    %fprintf('\tTPR: m=%d n=%d p=%d h=%d: %d\n',m,n,p,h,check);
                    assert(check<10^-T,'ERROR: does not pass inverse test check=%f',check);
                catch xy
                    fprintf('TPR: m=%d n=%d p=%d h=%d\n',m,n,p,h);
                    fprintf('TPR: a1=%d a2=%d b1=%d b2=%d\n',sum(a1),sum(a2),sum(b1),sum(b2));
                    disp(xy.message);
                    whos a1 a2 A b1 b2 B
                    %disp(full(A))
                    %disp(full(B))
                    %disp(full(B*A))
              end
	    catch AB
                fprintf('TPR: m=%d n=%d p=%d h=%d\n',m,n,p,h);
	    	disp(AB.message);
	    end

        end
    end
end

disp('end test');
