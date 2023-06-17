function add_var_shade_1tp(pt, r, col)
if numel(r) == 1; r = repmat(r,length(pt),1);end
plot_inputs = {'FaceAlpha',0.4,'EdgeColor' ,'none'};
if length(pt) == 3
    [X,Y,Z] = sphere;
    surf(r(1)*X + pt(1),r(2)*Y + pt(2),r(3)*Z + pt(3),repmat(permute(col,[3 1 2]),size(X,1),size(X,2),1), plot_inputs{:})
else
    theta = (0:0.1:2*pi)';
    patch(pt(1) + r(1).*cos(theta),pt(2) + r(2).*sin(theta),col, plot_inputs{:});
end
end