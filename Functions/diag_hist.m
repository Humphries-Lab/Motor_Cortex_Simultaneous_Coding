function diag_hist(X,max_val)

a=max_val*sqrt(2)/2;
b=a/sqrt(2)+0.005;

% rotation matrix
phi = -pi/4;
R = [cos(phi),-sin(phi);sin(phi),cos(phi)];

h=histogram(X,-a:0.001:a,'Normalization','probability','Visible','off');

% rotate and traslate datapoints
scaley=0.01;
tmpy=reshape([h.Values*scaley;h.Values*scaley],[],1);
tmpx=reshape([h.BinEdges;h.BinEdges],[],1);
coords3=[tmpx';[0 tmpy' 0]];
coords3 = R * coords3;

plot(coords3(1,:)+b,coords3(2,:)+b,'Color','c')

zer=[h.BinEdges;zeros(1,numel(h.BinEdges))];
coords2 = R * zer;
plot(coords2(1,:)+b,coords2(2,:)+b,'k')

scale=R * [0.001 00.001; 0 0.1*scaley];
plot(scale(1,:)+b,scale(2,:)+b,'r')

end