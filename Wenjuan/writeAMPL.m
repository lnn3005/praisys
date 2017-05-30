function writeAMPL(fileID,datatype,dataname,data)

% INPUT: fileID ~ file to be written
%        datatype ~ type of data to be written.
%                   1 = set
%                   11 = set of elements that contain 2 indices.
%                   2 = parameter. Must be number, array, or matrix
%        dataname ~ name of set or parameter
%        data ~ actual data to be printed
% OUTPUT: write data in an AMPL data file


% Set 
if datatype == 1
    printString = ['set ',dataname,' := ',num2str(data),';'];
    fprintf(fileID,'%s\n',printString);
end

% Set with the form (a,b) (c,d) (e,f) etc.
if datatype == 11
    data11 = ' ';
    for i = 1:size(data,1)
        data11 = strcat(data11,' (',num2str(data(i,1)),',',num2str(data(i,2)),')');
    end
    printString = ['set ',dataname,' := ',data11,';'];
    fprintf(fileID,'%s\n',printString);
end

% Parameter with 1 or 2 dimensions
if datatype == 2
    % If the parameter is a number
    if size(data,1) == 1 && size(data,2) == 1
        printString = ['param ',dataname,' := ',num2str(data),';'];
        fprintf(fileID,'%s\n',printString);        
    else % If the parameter is a matrix or an array    
        printString = ['param ',dataname,' : ',num2str(1:size(data,2)),':='];
        fprintf(fileID,'%s\n',printString);
        
        % Add row number to the parameter matrix
        newdata = [[1:size(data,1)]',data];
        fprintf(fileID,[repmat('%d\t', 1, size(newdata, 2)) '\n'], newdata');
        fprintf(fileID,'; \n');
    end   
end

% Parameter with 1 dimension
if datatype == 21 
    printString = ['param ',dataname,' := ',num2str(data),';'];
    fprintf(fileID,'%s\n',printString);            
end

% Parameter with 3 dimensions. Print along 3rd dimension 
if datatype == 23
    printString = ['param ',dataname,':='];
    fprintf(fileID,'%s\n',printString);
    
    for i=1:size(data,3)
        printString = ['[*,*,',num2str(i),'] : ',num2str(1:size(data,2)),' :='];
        fprintf(fileID,'%s\n',printString);
        newdata = [[1:size(data,1)]',data(:,:,i)];
        fprintf(fileID,[repmat('%d\t', 1, size(newdata, 2)) '\n'], newdata');
    end
    fprintf(fileID,'; \n');
end

fprintf(fileID,'\n');

end