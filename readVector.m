
function v = readVector(filename)
    % Read data from file
    data = load(filename);
    
    % Convert 0-based index to 1-based index (this assumes the indexing is consistent and starts at 0)
    indices = data(:, 1) + 1;
    
    % Extract vector values
    values = data(:, 2);
    
    % Initialize a zero vector of the appropriate length
    v = zeros(max(indices), 1);
    
    % Assign values based on indices
    v(indices) = values;
end