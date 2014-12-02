function whichb = uncatIndexArray(N)
%makes an index array for un-concatenation
%INPUT
%N:         1 x nb, dim of each block in concatenated array, e.g. A = [B1 B2 B3], N = [size(B1,2) size(B2,2) size(B3,2)]
%OUTPUT
%whichb:    1 x sum(N), indicates which block each row/col belongs to
%           use to index the concatenated array, e.g. B3 = A(:,whichb==3)

whichb = zeros(1,sum(N)); %which block
i = 0;
for bno = 1:length(N)
    I = i + (1:N(bno));
    whichb(I) = bno;
    i = I(end);
end