classdef opjRead
    %opjRead Read a .js file  
    % A method allowing users of Matlab to read in JavaSeis (.js) files
    % storing seismic data. This will return the file in a Matlab array and
    % give access to the grid information: eg, dimensions, offsets, and
    % intervals.
    
    %Note: This is wrong; need to actually read the .js file with javaseis
    %then convert to matlab, unfortunately. 
    
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
            function op = opjRead(filename)
                javaaddpath([pwd, '/jstest.jar'])
                javaaddpath([pwd, '/beta_javaseis.jar'])
                javaaddpath([pwd,'/edu_mines_jtk.jar'])
                import('beta.javaseis.array.ElementType.*');
                import('beta.javaseis.array.MultiArray.*');
                import('beta.javaseis.array.MultiArraySampleIterator.*');
                import('beta.javaseis.array.MultiArrayTraceIterator.*');
                import('beta.javaseis.array.TransposeType.*');
                import('beta.javaseis.io.Seisio');
                import ('beta.javaseis.grid.GridDefinition.*');
                import ('edu.mines.jtk.util.ParameterSet.*');
                import ('java.lang.System.*');
                
                % Create seisio
                tio = beta.javaseis.io.Seisio(filename);
                tio.open('r');
                shape=tio.getGridDefinition().getAxisLengths();
                dimensions = tio.getGridDefinition().getNumDimensions(); 
                elements=1;
                for i=1:dimensions
                    elements=elements*single(shape(i));
                end
                fid=fopen([filename,'/TraceFile']);
                op.matrix=fread(fid, elements, 'float');
                % reshape for different sizes:
                if tio.getGridDefinition().getNumDimensions() == 1
                    op.matrix=builtin('reshape',op.matrix,shape(1));
                elseif tio.getGridDefinition().getNumDimensions() == 2
                    op.matrix=builtin('reshape',op.matrix,shape(1), shape(2));
                elseif tio.getGridDefinition().getNumDimensions() == 3
                    op.matrix=builtin('reshape',op.matrix,shape(1), shape(2), shape(3));
                elseif tio.getGridDefinition().getNumDimensions() == 4
                    op.matrix=builtin('reshape',op.matrix,shape(1), shape(2), shape(3), shape(4));
                elseif tio.getGridDefinition().getNumDimensions() == 5
                    op.matrix=builtin('reshape',op.matrix,shape(1), shape(2), shape(3), shape(4), shape(5));
                end %reshaping
            end % function opjRead
        end % Methods
end % Classdef



% For an example, try: 
% opj= opjRead('/users/slic/troberso/demo_CP_code/test.js')