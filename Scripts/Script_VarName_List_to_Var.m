% Convert VarName_List into variable with names VarName_List
% Dynamically allocate the name of the variables
assert(all(size(VarName_List) == size(Var_List)));
for k = 1:length(VarName_List(:))
    VarName = VarName_List{k};   % Name of the variables
    VarVal = Var_List{k}.Value;   % Values of the variables
    eval([VarName,'=VarName_Range{k}{VarVal};']);
end
