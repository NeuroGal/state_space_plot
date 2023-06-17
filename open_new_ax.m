function new_ax = open_new_ax(ax_handle, pos_ax_units, axSize)
% new ax in according to pos in ax_handle, axSize is normalized units
% see this for the 3d version: https://www.mathworks.com/matlabcentral/answers/248362-screen-2d-projection-of-3d-plot
xl = get(ax_handle, 'XLim'); yl = get(ax_handle, 'YLim');
if length(pos_ax_units) == 3
    v = ax_handle.View; zl = get(ax_handle, 'ZLim');
    mat = viewmtx(v(1),v(2)) * makehgtform('scale',1./ax_handle.DataAspectRatio);
    edge_dots = [[xl;yl;zl],[xl(1);yl(2);zl(1)],[xl(2);yl(1);zl(1)]];
    pts_new = mat * [[pos_ax_units',edge_dots];ones(1,5)];
    pts_new_2d = pts_new(1:2,:)./pts_new(4,:);
    pos_ax_units_2d = pts_new_2d(:,1)';
    xl = pts_new_2d(1,4:5); yl = pts_new_2d(2,2:3);
else
    pos_ax_units_2d = pos_ax_units;
end
xycNorm = (pos_ax_units_2d - [xl(1),yl(1)])./[range(xl),range(yl)]; %normalized to axis
set(ax_handle, 'Units', 'Normalize'); axpos = get(ax_handle,'position');
xycFigNorm = axpos(1:2) + axpos(3:4).*xycNorm; %normalized to figure 
new_pos = [xycFigNorm, axSize];
new_ax = axes('Units','Normalize','Position',new_pos, 'color', 'none');
end