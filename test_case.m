%fid=fopen("matrices/bcsstk03.mtx");
%e=textscan(fid,'%f %f %f','CommentStyle','%');
%fclose(fid);
%out = cell2mat(e);

%n = out(1, 1); % Number of rows
%m = out(1, 2); % Number of columns
%nnz_vals = out(1, 3); % Number of non-zero elements
%row_indices = out(2:end, 1);
%col_indices = out(2:end, 2);
%values = out(2:end, 3);
%toc

%lower_mat = sparse(out(2:end, 1), out(2:end, 2), out(2:end, 3), out(1, 1), out(1, 2), out(1, 3));

%sparse_mat = lower_mat + lower_mat' - diag(diag(lower_mat));
%tic
% Reorder
%display("Reorder")
%perm = symamd(sparse_mat); % Find permutation for sparsity
%reordered_mat = sparse_mat(perm, perm); % Reorder the matrix
lower_mat = mmread("matrices/bcsstk03.mtx");
tic

display("Calculate")
mat = chol(lower_mat, 'lower');

display(nnz(mat))
toc

mmwrite("result_mat.mtx", mat, '', 'real', 6);

%fid=fopen("result.mtx");
%e=textscan(fid,'%f %f %f','CommentStyle','%');
%out_ = cell2mat(e);
%res = sparse(out_(2:end, 1), out_(2:end, 2), out_(2:end, 3), out_(1, 1), out_(1, 2), out_(1, 3));


%testmat = mat-res;
%display(norm(testmat, Inf)/norm(mat, Inf));
%display(nnz(res));
%display(nnz(mat));

%tic
%fid=fopen("result.mtx");
%e=textscan(fid,'%f %f %f','CommentStyle','%');
%fclose(fid);
%out = cell2mat(e);
%lower_mat = sparse(out(2:end, 1), out(2:end, 2), out(2:end, 3), out(1, 1), out(1, 2), out(1, 3));

%display(nnz(lower_mat));
%display(nnz(mat));

%display(max(abs(lower_mat - mat),[],"all","linear"));
%display(max(abs(lower_mat - mat),[],"all")/mean(mat,"all"));
%errors =abs(lower_mat - mat)/max(lower_mat,mat);
%errors =abs(lower_mat - mat)/mat;
%subplot(1, 2, 1), spy(mat), title('Program');
%subplot(1, 2, 2), spy(lower_mat), title('Mine');
%lower_mat = lower_mat(1:100000);
%mat = mat(1:100000);
%pause();

%disp(max(errors, [], "all"));
%toc
