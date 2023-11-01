% display network parameters in an image (2D)
% with one axes as neuron input, and other as output
% p is the parameter to be displayed
function P = display_network_param(X,Xn,p,N,fig)
% create sparse matrix from matlab variables
if(nargin>1)
    P = sparse(X,Xn,p,N,N); % i->X, j-> Xn, val ->p
else
    P = X;
end

if(nargin>=5)
    if(isempty(fig))
        figure(3);
    else
    if(isgraphics(fig,'figure'))
    figure(fig);
    else(isgraphics(fig,'subplot'))
        subplot(fig);
    end
    end
end
green_fraction = full(max(P(:))/(max(P(:))-min(P(:))));
green_id = uint16(256*green_fraction);
greenColorMap = [zeros(1, 256-green_id), linspace(0, 1,green_id)];
redColorMap = [linspace(1, 0, 256-green_id), zeros(1, green_id)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
colormap(colorMap);
imagesc(P);colorbar;
xlabel('Out to Neurons');
ylabel('Input');
