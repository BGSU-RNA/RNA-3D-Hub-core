function checkMatchingMatrix(M)

    valid = 1;

    [x,y] = size(M);
    if x ~= y
        disp('Matrix is not square');
        valid = 0;
    end
    
    if ~isequal(triu(M), tril(M)') || ~isequal(tril(M), triu(M)')
        disp('Matrix is not symmetric');
        valid = 0;        
    end

    if ~isempty(find(diag(M) ~= 0, 1))
        disp('Non-zero entries on the diagonal');
        valid = 0;        
    end

    if ~isempty(find(M < 0, 1))
        disp('Negative values found');
        valid = 0;        
    end    

    if ~isempty(find(diag(M) ~= 0, 1))
        disp('Non-zero entries on the diagonal');
        valid = 0;        
    end        
    
    if valid == 1
        disp('Matrix is square and symmetric');
    end
            
end