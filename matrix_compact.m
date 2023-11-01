%from
%https://in.mathworks.com/matlabcentral/answers/127013-how-can-i-eliminate-unwanted-zeros-from-a-matrix 
function B =matrix_compact(A)
 T = A~=0;
   n = sum(T,2);
   m = max(n);
   B = nan(size(A,1),m);
   for k = 1:size(A,1)
     B(k,m-n(k)+1:m) = A(k,T(k,:)); % <-- B is the result
   end