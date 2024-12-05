tic
fid=fopen("ct20stif.mtx");
e=textscan(fid,'%f %f %f','CommentStyle','%');
out = cell2mat(e);
mat = chol(sparse(out(2:end, 1), out(2:end, 2), out(2:end, 3), out(1, 1), out(1, 2), out(1, 3)), 'lower');
toc

fid=fopen("result.mtx");
e=textscan(fid,'%f %f %f','CommentStyle','%');
out_ = cell2mat(e);
res = sparse(out_(2:end, 1), out_(2:end, 2), out_(2:end, 3), out_(1, 1), out_(1, 2), out_(1, 3));


testmat = mat-res;
display(norm(testmat, Inf)/norm(mat, Inf));
display(nnz(res));
display(nnz(mat));
