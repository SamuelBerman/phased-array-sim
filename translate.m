function [X, Y, Z] = translate(Xi, Yi, Zi, tx, ty, tz)
    h=numel(Xi)/2;
    X=Xi;
    Y=Yi;
    Z=Zi;
    X(1:h)=Xi(1:h)+tx;
    Y(1:h)=Yi(1:h)+ty;
    Z(1:h)=Zi(1:h)+tz;
end