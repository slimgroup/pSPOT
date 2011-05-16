function distHeaderWrite(file_name, n, d, o, l, u)
%distHeaderWrite writes header to xml file



% Check number of arguments
assert(nargin >= 2,'There must be at least 2 input arguments');

dims = length(n);

% Set the default values in case we have less tha 6 arguments
if(nargin < 6)    
    for i=1:dims
    u(i) = {strcat('u',int2str(i))};
    end
    if(nargin < 5)
        for i=1:dims
        l(i) = {strcat('l',int2str(i))};
        end
    end
    if(nargin < 4)
        for i=1:dims
        o(i) = 0;
        end
    end
    if(nargin < 3)
        for i=1:dims
        d(i) = 1;
        end
    end
end % if

docNode = com.mathworks.xml.XMLUtils.createDocument... 
    ('header')

docRootNode = docNode.getDocumentElement;

% Writing n to xml
thisElement = docNode.createElement('n');
thisElement.appendChild... 
    (docNode.createTextNode(sprintf('%i',n(1))));
for i=2:dims
    thisElement.appendChild... 
    (docNode.createTextNode(sprintf(' %i',n(i))));
end % for
docRootNode.appendChild(thisElement);

% Writing d to xml
thisElement = docNode.createElement('d'); 
thisElement.appendChild... 
    (docNode.createTextNode(sprintf('%i',d(1))));
for i=2:dims
    thisElement.appendChild... 
    (docNode.createTextNode(sprintf(' %i',d(i))));
end % for
docRootNode.appendChild(thisElement);

% Writing o to xml
thisElement = docNode.createElement('o'); 
thisElement.appendChild... 
    (docNode.createTextNode(sprintf('%i',o(1))));
for i=2:dims
    thisElement.appendChild... 
    (docNode.createTextNode(sprintf(' %i',o(i))));
end % for
docRootNode.appendChild(thisElement);

% Writing l to xml
thisElement = docNode.createElement('l'); 
thisElement.appendChild... 
    (docNode.createTextNode(sprintf('%s',l{1})));
for i=2:dims     
    thisElement.appendChild... 
    (docNode.createTextNode(sprintf(' %s',l{i})));
end % for
docRootNode.appendChild(thisElement);

% Writing u to xml
thisElement = docNode.createElement('u');
thisElement.appendChild... 
    (docNode.createTextNode(sprintf('%s',u{1})));
for i=2:dims    
    thisElement.appendChild... 
    (docNode.createTextNode(sprintf(' %s',u{i})));
end % for
docRootNode.appendChild(thisElement);

% Setting the xml file name
xmlFileName = [file_name,'.xml'];
xmlwrite(xmlFileName,docNode);
type(xmlFileName);

end % function