% testing window functions
M=[12 14];
P=4;
H=3;

disp('start test');

for m=M(1):M(2)
    for p=2:P
        for h=0:H
    
            [A a1 a2]=pWindow.fd_window_1Dfor(m,p,h);
            [B b1 b2]=pWindow.fd_window_1Dbck(m,p,h);

            try
                n=m+(p-1)*2*h;
                x0=rand(m,1);
                y=A*x0;
                x1=B*y;
                check=sum(abs((x1-x0)'));
                fprintf('\tm=%d n=%d p=%d h=%d: %d\n',m,n,p,h,check);
            catch me
                fprintf('m=%d n=%d p=%d h=%d\n',m,n,p,h);
                fprintf('a1=%d a2=%d b1=%d b2=%d\n',sum(a1),sum(a2),sum(b1),sum(b2));
                whos a1 a2 A b1 b2 B
            end


        end
    end
end

disp('end test');
