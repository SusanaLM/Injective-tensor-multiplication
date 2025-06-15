function hyperdet = hyperdeterminant(T)

if ~isequal(size(T), [2 2 2])
    error('Input tensor must be of size 2x2x2.');
end

% Extract components
t000 = T(1,1,1);
t001 = T(1,1,2);
t010 = T(1,2,1);
t011 = T(1,2,2);
t100 = T(2,1,1);
t101 = T(2,1,2);
t110 = T(2,2,1);
t111 = T(2,2,2);

% Compute hyperdeterminant according to Cayley's formula
hyperdet = t000^2 * t111^2 + t001^2 * t110^2 + t010^2 * t101^2 + t100^2 * t011^2 ...
    - 2 * (t000*t001*t110*t111 + t000*t010*t101*t111 + t000*t011*t101*t110 + ...
           t001*t010*t101*t110 + t001*t011*t100*t111 + t010*t011*t100*t101) ...
    + 4 * (t000*t011*t101*t110 + t001*t010*t100*t111);
end