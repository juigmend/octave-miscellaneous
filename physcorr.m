function C = physcorr(a,b)
%PHYSCORR Physical Correlation
%
%   PHYSCORR(a,b) where a and b are row vectors of equal length.
%   Each of these vectors contain a curve.
%
%   Physical Correlation does not remove the mean and thus is useful,
%   for example, for signals that are equally offset
%   (c.f.Bruderer, McKinney & Kohlrausch, 2012).
%   
%   Juan Ignacio Mendoza - 2016

C = (a * b')/ sqrt( sum(a.^2) * sum(b.^2) );

end
