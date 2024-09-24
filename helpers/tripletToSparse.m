function S = tripletToSparse(filename)
    % Read triplet data from file
    data = load(filename);
    
    % Extract row, column, and value information
    rows = data(:, 1) + 1; % Convert to 1-based indexing
    cols = data(:, 2) + 1; % Convert to 1-based indexing
    values = data(:, 3);
    
    % Find matrix dimensions
    m = max(rows);
    n = max(cols);
    
    % Convert triplet form to full matrix
    S = sparse(rows, cols, values, m, n);
    
    % Mirror upper-right into lower-left for symmetric matrix
    S = S + S' - diag(diag(S));
end

