A = {};
for k=1:10
%     A{k} = randn(2^20/(2^3), 3653);
    A{k} = randn(2^20/(2^8), 3653);
end
s = whos('A');
A_size = s.bytes/(1024^3);  % 1024^3: in GBytes