function ymin = y_min_stat(geo)

[ay, cy, cx, L, l, A, B, C] = unpack_geo(geo);

ymin_t = max((-A(2)/B(2)),(-A(3)/B(3)));

if (L - ymin_t)< 0
    ymin = NaN;
else
    ymin = ymin_t;
end
