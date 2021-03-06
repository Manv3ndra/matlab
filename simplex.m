clc
clear

fprintf("Code run\n");
constants = {'x1','x2','x3','x4','s1','s2','s3','sol'};
%constraints
a = [1 4 8 6; 4 1 2 1; 2 3 1 2];
%constraint solution
b = [11;7;2];
%objective function
C = [4 6 3 1];
%objective function size
[m,n] = size(a);
%initialize s variable
s = eye(m);
%main tableau
A = [a s b];
%calculate initial cost
cost = zeros(1,n+m+1);
cost(1:n) = C;
%basic variables
bv = n+1:1:n+m;
%zj - cj
zjcj = cost(bv)*A - cost;
%main tableau
zcj = [zjcj ; A];
smpTb = array2table(zcj);
smpTb.Properties.VariableNames(1:n+m+1) = constants;
disp(smpTb);

while true
    %check if zj - cj > 0
    if(any(zjcj < 0))
        fprintf('The current BFS is not optimal\n');
        %zj - cj without solution value
        zc = zjcj(1:n+m);
        %get entering variable and its column index
        [entering_val,pvt_col] = min(zc);
        if(all(A(:, pvt_col) < 0))
            fprintf("LPP is unbounded\n");
        else
            %solution column
            sol(:) = A(:,n+m+1);
            %entering column
            column = A(:,pvt_col);
            %calculate ratios
            ratio = zeros(m,1);
            for i =1:m
                if(column(i)>0)
                    ratio(i) = sol(i)./column(i);
                else
                    ratio(i) = inf;
                end
            end
            %get leaving variable and its row index
            [leaving_val, pvt_row] = min(ratio);
            %update basic variables
            bv(pvt_row) = pvt_col;
            %get pivot key
            pvt_key = A(pvt_row,pvt_col);
            %row operation -> divide row by pivot key
            A(pvt_row,:) = A(pvt_row,:)./pvt_key;
            %convert all values of pivot column to 0 except key with row operation R[i] = R[i] - R[i,pvt_col]*R[pivot_row]
            for i = 1:m
                if i ~= pvt_row
                    A(i,:) = A(i,:) - A(i,pvt_col).*A(pvt_row,:);
                end
            end
            %calculate zj - cj
            zjcj = zjcj - zjcj(pvt_col).*A(pvt_row,:);
            zcj = [zjcj; A];
            table = array2table(zcj);
            table.Properties.VariableNames(1:n+m+1) = constants;
            disp(table);
        end
    else
        fprintf('The current BFS is optimal\n');
        break;
    end
end

%print solution
bv = [n+m+1, bv];
for i=1:m+1
    var = bv(i);
    fprintf('%s = %f\n',constants{var},zcj(i,n+m+1));
end
