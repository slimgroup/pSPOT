classdef opjJavaToMatlab
    %opjRead Read a .js file  
    % A method allowing users of Matlab to read in JavaSeis (.js) files
    % storing seismic data. This will return the file in a Matlab array and
    % give access to the grid information: eg, dimensions, offsets, and
    % intervals.
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        properties( SetAccess = protected )
           matrix
        end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Methods - Public
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        methods   
            %Constructor
            function op = opjJavaToMatlab(javaArray)
                import('beta.javaseis.array.MultiArray.*');
                length = javaArray.getTotalElementCount();
                op.matrix = zeros(length, 1);
                iter = beta.javaseis.array.MultiArraySampleIterator( javaArray );
                i=1;
                while iter.hasNext()
                    op.matrix(i) = iter.getNextFloat();%does JavaSeis write in floats?I suppose we should check type.
                    i=i+1;
                end
                op.matrix=reshape(op.matrix, 20,5);             
            end % function opjRead
        end % Methods
end % Classdef
