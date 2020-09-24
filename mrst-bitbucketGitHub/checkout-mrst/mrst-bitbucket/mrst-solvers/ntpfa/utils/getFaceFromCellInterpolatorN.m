function [P, faceNo, cellNo, w] = getFaceFromCellInterpolatorN(G, state, rock, ndof)
    if nargin == 2
        ndof = G.cells.num;
    end
    
    if ~isfield(state, 'N_f')
        [P, faceNo, cellNo, w] = getFaceFromCellInterpolator(G, rock, ndof);
        return
    end
    % Get nodes
    faceNo = G.cells.faces(:, 1);
    pos = G.cells.facePos;
    cellNo = rldecode((1:G.cells.num)', diff(pos));
    
    L = (G.cells.centroids(cellNo, :) - G.faces.centroids(faceNo, :));
    % One over distance weighting
    w = 1./sqrt(sum(L.^2, 2));
    
    N_f = state.N_f;
    N_c = state.N_c;
    
    N = G.faces.neighbors;
    intx = all(N > 0, 2);
    [bfaces, bcells] = boundaryFaces(G);

%     p_l = zeros(G.faces.num, 1);
%     p_r = zeros(G.faces.num, 1);
%     p_l(intx) = state.pressure(N(intx, 1));
%     p_r(intx) = state.pressure(N(intx, 2));
%     p_l(~intx) = state.pressure(bcells);
%     l_r(~intx) = state.boundaryPressure;
    
    N_f_t = accumarray(faceNo, N_f);
    w = N_c./N_f_t(faceNo);

%     isBF = false(G.faces.num, 1);
%     isBF(bfaces) = true;
%     w = ;
%     w = N_c
    
    % Divide by sum of weights for partition of unity
%     sumw = accumarray(faceNo, w);
%     w = w./sumw(faceNo);
    
    % n_n by n_c matrix constructing values on nodes by cells
    P = sparse(faceNo, cellNo, w, G.faces.num, ndof);
end

function K = getK(rock, cell, dim)
    k = rock.perm(cell, :);
    switch numel(k)
        case 1
            % Scalar perm
            K = k;
        case dim
            % Diagonal tensor
            K = diag(k);
        case 3*(dim - 1)
            % Full symmetric tensor
            if dim == 2
                K = [k(1), k(2); ...
                     k(2), k(3)];
            else
                K = [k(1), k(2), k(3); ...
                     k(2), k(4), k(5); ...
                     k(3), k(5), k(6)];
            end
        otherwise
            error('What sorcery is this?!');
    end
end
