%% Get euclidian distance between vector matrix
function D = euclidDistance(R)
[~,Ndim] = size(R);
D = (R(:,1)-R(:,1)').^2;
    for i = 2:Ndim
        D = D + (R(:,i)-R(:,i)').^2;
    end
end