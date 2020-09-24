function op = getInterpolationOperator(Pcell, Pface, Pnode, globalIndices, types, N)

%     Pcell = speye(G.cells.num);
    isCell = types == 1;
    isFace = types == 2;
    isNode = ~isCell & ~isFace;
    
    P = sparse(numel(globalIndices), size(Pcell, 2));
    if any(isCell)
        P(isCell,:) = Pcell(globalIndices(isCell), :);
    end
    if any(isFace)
        P(isFace,:) = Pface(globalIndices(isFace), :);
    end
    if any(isNode)
        P(isNode,:) = Pnode(globalIndices(isNode), :);
    end
    [ii, jj, vv] = find(P);
    keep = N(ii, 1) == jj;
    if size(N, 2) > 1
        keep = keep | N(ii, 2) == jj;
    end
    P_active = sparse(ii(keep), jj(keep), vv(keep), size(P, 1), size(P, 2));
    op = struct('P', P, 'P_active', P_active, 'P_passive', P - P_active);
    
 
%     op = @(pc) operator(pcell, pface, pnode, isCell, isFace, isNode, pc);
end

function val = operator(P_cell, P_face, P_node, isCell, isFace, isNode, v_cell)
    val = zeros(size(isCell));
    if isa(v_cell, 'ADI')
        val = double2ADI(val, v_cell);
    end
    if any(isCell)
        val(isCell) = P_cell*v_cell;
    end
    
    if any(isFace)
        val(isFace) = P_face*v_cell;
    end
    
    if any(isNode)
        val(isNode) = P_node*v_cell;
    end
end