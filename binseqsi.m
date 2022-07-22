function [S, d, c, f, l, m] = binseqsi(a,b,L,r)
%BINSEQSI Binary Sequences Similarity version 2
%
%   BINSEQSI(a,b,L,r) where a and b are row vectors of equal or different
%   length. These vectors contain indexes of either value in a binary sequence,
%   therefore they should not have repeated elements within each one.
%
%   L is the length of both sequences.
%   r = 0 turns display of results off.
%   r = 1 turns display of results on (default).
%
%   [S, d, c, f, l, m] = binseqsi(a,b,L,r)
%
%       S: Similarity value between 0 and 1;
%          The greater the value, the more similar the sequences are.
%
%       d: Distance of paired elements.
%
%       c: Closeness of paired elements.
%
%       f: Fraction of paired elements.
%
%       l: Lag of paired elements.
%
%       m: Cell array containing the paired elements of a and b.
%
%   Juan Ignacio Mendoza - 2016

timervalue = tic;

% check input arguments:

err_r = ('binseqsi ERROR: r should be 1 or 0');

if r > 1
    disp(err_r)
    return
elseif r < 0
    disp(err_r)
    return
end

% check that the vectors do not contain zeroes:

if isempty(find([a, b] == 0)) == 0
    disp 'binseqsi ERROR: vectors should not contain zeroes.';
    return
end

% check that L is not smaller than the larger of the indexes:

a = sort(a);
b = sort(b);

largerindex = max(a(end),b(end));

if L < largerindex
    disp 'binseqsi ERROR: Length of sequences should not be smaller than the greatest element of both index vectors.';
    return
end

size_a = length(a);
size_b = length(b);
maxsize = max(size_a, size_b);
minsize = min(size_a, size_b);

% check that the vectors do not contain duplicates:

errdups = ('binseqsi ERROR: vectors should not contain duplicates.');

if length(unique(a)) ~=  size_a
    disp(errdups);
    return
end

if length(unique(b)) ~=  size_b
    disp(errdups);
    return
end

% make reference matrices for each vector:

id_a = repmat(a',1,size_b);
id_b = repmat(b,size_a,1);

% =========================================================================
% Pair each element of vector b with an element of vector a
% From the paired elements compute Lag, Closeness, Fraction of
% Paired Elements and Similarity.

% distance matrix:

for i_1 = 1:size_a %
    for i_2 = 1:size_b
        distmat(i_1,i_2) = -(a(i_1)-b(i_2));
    end
end

absdistmat = abs(distmat);

% Make a logical matrix indicating the minimum value(s) of each column:

for i_1 = 1:size_a
    mincols(i_1,:) = absdistmat(i_1,:) == min(absdistmat);
end

% Make a logical matrix indicating the minimum value(s) of each row:

for i_1 = 1:size_b
    minrows(:,i_1) = absdistmat(:,i_1) == min(absdistmat,[],2);
end

allmins = mincols .* minrows; % intersection
allmins(allmins == 0) = NaN;

% find to which elements of original a and b do these minima
% correspond:

from_a = allmins.*id_a;
m{1,1} =  (from_a(isfinite(from_a)))';

from_b = allmins.*id_b;
m{2,1} = (from_b(isfinite(from_b)))';

% compute measures:

l = (nanmean(nanmean(allmins(:) .* distmat(:)))); % lag
d = (nanmean(nanmean(allmins(:) .* absdistmat(:)))); % distance
minunel = min(size(unique(m{1,1}),2),size(unique(m{2,1}),2));

c = (1-(d/L)); % paired elements closeness
f = minunel/maxsize; % fraction of paired elements
S = c*f;

% =========================================================================
% Display results

if r == 1
    
    toc(timervalue)
    disp --------------------------------------------------
    disp ' BINARY SEQUENCES SIMILARITY '
    disp ' '
    fprintf(' distance = %f \n', d);
    fprintf(' closeness = %f \n', c);
    fprintf(' lag = %f \n', l);
    fprintf(' fraction of paired elements = %f \n', f);
    disp ' '
    fprintf(' Similarity = %f \n', S);
    disp ' '
    disp 'Paired elements:'
    disp ' '
    disp(m{1,1})
    disp(m{2,1})
    disp --------------------------------------------------
    
elseif r == 0
    return
end

end