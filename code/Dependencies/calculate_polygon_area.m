
function area = calculate_polygon_area(vertices)
    n = length(vertices);
    area = 0;
    for i = 1:n
        j = mod(i, n) + 1;
        area = area + (vertices(i, 1) * vertices(j, 2) - vertices(j, 1) * vertices(i, 2));
    end
    area = abs(area) / 2;
end