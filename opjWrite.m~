classdef opjWrite
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
            function op = opjWrite(filename, array)
                javaaddpath([pwd, '/beta_javaseis.jar'])
                javaaddpath([pwd,'/edu_mines_jtk.jar'])
                import('beta.javaseis.array.ElementType.*');
                import('beta.javaseis.array.MultiArray.*');
                import('beta.javaseis.array.MultiArraySampleIterator.*');
                import('beta.javaseis.array.MultiArrayTraceIterator.*');
                import('beta.javaseis.array.TransposeType.*');
                import('beta.javaseis.io.Seisio.*');
                import ('beta.javaseis.grid.GridDefinition.*');
                import ('edu.mines.jtk.util.ParameterSet.*');
                import ('java.lang.System.*');
                

                grid = beta.javaseis.grid.GridDefinition.getDefault( ndims(array), size(array) );
                sio = beta.javaseis.io.Seisio( filename, grid );
                sio.create();
                
                sa = beta.javaseis.array.MultiArray(ndims(array), beta.javaseis.array.ElementType.DOUBLE, size(array));
                iter = beta.javaseis.array.MultiArraySampleIterator(sa);
                array = builtin( 'reshape', array, numel(array), 1 );
                index=1;
                while iter.hasNext(),
                    iter.next();
                    iter.putDouble(array(index)); %Error here - null pointer, and I don't know why :(
                    index = index + 1;
                end
                
                sio.writeMultiArray( sa, zeros(1 , ndims(array)) ); % multiarray must have at least 2 and no more than 5 dimensions
                sio.close();   
                
            end % function opjRead
        end % Methods
end % Classdef