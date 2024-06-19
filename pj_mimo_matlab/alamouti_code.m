function codedData = alamouti_code(data)
    if mod(length(data), 2) ~= 0
        error('The length of the input data must be even for Alamouti coding.');
    end
    data = reshape(data, 2, []);
    codedData = zeros(2, size(data, 2) * 2);
    for i = 1:size(data, 2)
        s1 = data(1, i);
        s2 = data(2, i);
        codedData(:, 2*i-1) = [s1; s2];        % first time slot
        codedData(:, 2*i) = [-conj(s2); conj(s1)];  % second time slot
    end
end
