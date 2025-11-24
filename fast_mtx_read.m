function A = fast_mtx_read(filename)
    % Fastest pure-Octave coordinate .mtx reader

    fid = fopen(filename, "r");
    if fid < 0, error("Cannot open file"); end

    % Skip comments and read header
    line = fgetl(fid);
    while isempty(line) || line(1) == '%'
        line = fgetl(fid);
    end

    % Read size line
    sz = sscanf(line, "%d %d %d");
    rows = sz(1); cols = sz(2); nnz = sz(3);

    % Now read ALL remaining lines AS RAW TEXT in a single call
    raw = fread(fid, "*char")';
    fclose(fid);

    % Fastest parser: one vectorized sscanf
    data = sscanf(raw, "%d %d %f");

    % Reshape into columns: (i j v)
    data = reshape(data, 3, nnz);

    i = data(1, :)';
    j = data(2, :)';
    v = data(3, :)';

    % Construct sparse matrix
    A = sparse(i, j, v, rows, cols);
end