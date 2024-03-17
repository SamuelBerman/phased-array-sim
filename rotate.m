function [X, Y, Z] = rotate(Xi, Yi, Zi, tx, ty, tz)
    h=numel(Xi)/2;
    X=Xi;
    Y=Yi;
    Z=Zi;
    out = rotx(tx)*[Xi(1:h); Yi(1:h); Zi(1:h)];
    out = roty(ty)*out;
    out = rotz(tz)*out;
    X(1:h)=out(1,:);
    Y(1:h)=out(2,:);
    Z(1:h)=out(3,:);
end