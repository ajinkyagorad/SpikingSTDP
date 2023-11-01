% Check performance using gpu

x = gpuArray.rand(5000,5000);
tic
for i = 1:1000
    y = sum(x,1);
end
toc