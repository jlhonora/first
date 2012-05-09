function [y] = triang(x)

y = (x>-1).*(x<1).*(1-abs(x));