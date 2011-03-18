classdef opjMatlabToJava
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
            
            function op = opjMatlabToJava(matlabArray)
                sa = beta.javaseis.array.MultiArray.factory(ndims(matlabArray), ...
                    beta.javaseis.array.ElementType.DOUBLE, ...
                    numel(matlabArray), size(matlabArray));
                        %might need to adjust for different types
                iter = beta.javaseis.array.MultiArraySampleIterator(sa);
                matlabVector = builtin( 'reshape', matlabArray, numel(matlabArray), 1 );
                index=1;
                while iter.hasNext(),
                    iter.next();
                    iter.putDouble(matlabVector(index));
                    index = index + 1;
                end
                op.matrix=sa;
            end % function opjRead
        end % Methods
end % Classdef
