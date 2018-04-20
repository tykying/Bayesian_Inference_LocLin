function [SpDis_ranges_min, SpDis_ranges_max] = LocateBoundaries(x, y)
    SpDis_ranges_min = floor(min([x(:), y(:)]));
    SpDis_ranges_max = ceil(max([x(:), y(:)]));
end