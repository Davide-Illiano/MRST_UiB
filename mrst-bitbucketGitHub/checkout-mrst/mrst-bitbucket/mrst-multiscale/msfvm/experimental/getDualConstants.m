function DG = getDualConstants(CG, DG)
    G = CG.parent;
    is3d = CG.griddim == 3 && G.cartDims(3) > 1;
    [~, DG] = createPermutationMatrix(sparse(G.cells.num, G.cells.num), DG, CG, 'Speedup', is3d);

    P = double(DG.P);
    X = restrictOperator(CG, DG, DG.N);
    DG.X = X*P;
end

function X = restrictOperator(CG, DG, Nf)
    %Returns a matrix representing the restriction operator
    %Should be nxf big where n is coarse nodes and f total fine nodes
    %Permute the partition to the new index space
    permuted_partition = DG.P*CG.partition;
    xind = zeros(1,Nf);
    yind = zeros(1,Nf);
    pos = 1;
    for coarse = 1:CG.cells.num
       %find indices of corresponding coarse nodes
       indices = find(permuted_partition == coarse);
       %insert the correct cells into the arrays of indices
       M = length(indices);
       yind(pos:(pos + M-1)) = indices;
       xind(pos:(pos + M-1)) = coarse;
       pos = pos + M;
    end
    X = sparse(xind, yind, 1, Nf, Nf) > 0;
    return
end
