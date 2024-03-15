function [S, Tx, Ty, Tz] = dTrilinears3d(T, x, y, z, hx, hy, hz, bc)
    % function [P, Tx, Ty, Tz] = dTrilinears(T,x,y,z,hx,hy,hz,bc)
    % 
    % Uses trilinear particle in cell interpolation. Also returns
    % derivatives (hence the d in the name). 
    %
    %
    % T is the density. x is x+dx (final spatial x coordinate)
    %
    % hx, hy, hz are the dimensions of the cell along each coordinate (length, width, height). 
    %
    %
    % Use linear interpolation to evaluate (approximate) the values of T (density) in the
    % grid (x,y,z)

    % That = T(x,y,z) to be moved, necessary only for the derivative
    %
    % S is an operator that moves the mass according to the velocity field
    % and then redistributes the mass to cell centers (8 neighbors in 3D). 
    % Also returns derivatives called Tx Ty Tz which are used later in the
    % gradient descent for velociy. 
    %
    % Output: S, moves T, T^n+1 = S'*T^n;
    %         Tx,Ty,Tz - derivatives of (S'*T^n) w.r.t. x, y and z coordinate
    %% 
    %% 3D verison

    if nargin < 8
        %bc = 'closed';
        bc = 'open';
    end

    [m, n, s] = size(x);
    N1 = n * m;
    N = N1 * s;

    % Convert x and y to the coordinate system 1:m, 1:n, 1:s
    % x = 1 + x/h1 + 1/2; y = 1+ y/h2 + 1/2; z = 1+ z/h3 + 1/2;
    x = x ./ hx; y = y ./ hy; z = z ./ hz;
    x = x + 1/2; y = y + 1/2; z = z + 1/2;

    xd = floor(x(:)); xp = x(:) - xd;
    yd = floor(y(:)); yp = y(:) - yd;
    zd = floor(z(:)); zp = z(:) - zd;

    switch bc
        case 'closed'
            %CLOSED VERSION dtri1 (use this one)
            %
            ind1 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind2 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind3 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind4 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind5 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind6 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind7 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);
            ind8 = find(1 <= xd & xd < m & 1 <= yd & yd < n & 1 <= zd & zd < s);

        case 'open'
            % OPEN VERSION dtri3 (use this one)
            %
            ind1 = find(1 <= xd & xd < m + 1 & 1 <= yd & yd < n + 1 & 1 <= zd & zd < s + 1);
            ind2 = find(0 <= xd & xd < m & 1 <= yd & yd < n + 1 & 1 <= zd & zd < s + 1);
            ind3 = find(1 <= xd & xd < m + 1 & 0 <= yd & yd < n & 1 <= zd & zd < s + 1);
            ind4 = find(0 <= xd & xd < m & 0 <= yd & yd < n & 1 <= zd & zd < s + 1);
            ind5 = find(1 <= xd & xd < m + 1 & 1 <= yd & yd < n + 1 & 0 <= zd & zd < s);
            ind6 = find(0 <= xd & xd < m & 1 <= yd & yd < n + 1 & 0 <= zd & zd < s);
            ind7 = find(1 <= xd & xd < m + 1 & 0 <= yd & yd < n & 0 <= zd & zd < s);
            ind8 = find(0 <= xd & xd < m & 0 <= yd & yd < n & 0 <= zd & zd < s);

    end

    jki = xd(ind1) + (yd(ind1) - 1) * m + (zd(ind1) - 1) * N1;
    j1ki = xd(ind2) + 1 + (yd(ind2) - 1) * m + (zd(ind2) - 1) * N1;
    jk1i = xd(ind3) + yd(ind3) * m + (zd(ind3) - 1) * N1;
    j1k1i = xd(ind4) + 1 + yd(ind4) * m + (zd(ind4) - 1) * N1;

    jki1 = xd(ind5) + (yd(ind5) - 1) * m + zd(ind5) * N1;
    j1ki1 = xd(ind6) + 1 + (yd(ind6) - 1) * m + zd(ind6) * N1;
    jk1i1 = xd(ind7) + yd(ind7) * m + zd(ind7) * N1;
    j1k1i1 = xd(ind8) + 1 + yd(ind8) * m + zd(ind8) * N1;

    A1 = sparse(jki, ind1, (1 - xp(ind1)) .* (1 - yp(ind1)) .* (1 - zp(ind1)), N, N);
    A2 = sparse(j1ki, ind2, (xp(ind2)) .* (1 - yp(ind2)) .* (1 - zp(ind2)), N, N);
    A3 = sparse(jk1i, ind3, (1 - xp(ind3)) .* (yp(ind3)) .* (1 - zp(ind3)), N, N);
    A4 = sparse(j1k1i, ind4, (xp(ind4)) .* (yp(ind4)) .* (1 - zp(ind4)), N, N);
    A5 = sparse(jki1, ind5, (1 - xp(ind5)) .* (1 - yp(ind5)) .* (zp(ind5)), N, N);
    A6 = sparse(j1ki1, ind6, (xp(ind6)) .* (1 - yp(ind6)) .* (zp(ind6)), N, N);
    A7 = sparse(jk1i1, ind7, (1 - xp(ind7)) .* (yp(ind7)) .* (zp(ind7)), N, N);
    A8 = sparse(j1k1i1, ind8, (xp(ind8)) .* (yp(ind8)) .* (zp(ind8)), N, N);

    S = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8;

    %%
    if nargout > 1

        A1 = sparse(jki, ind1, T(ind1), N, N);
        A2 = sparse(j1ki, ind2, T(ind2), N, N);
        A3 = sparse(jk1i, ind3, T(ind3), N, N);
        A4 = sparse(j1k1i, ind4, T(ind4), N, N);
        A5 = sparse(jki1, ind5, T(ind5), N, N);
        A6 = sparse(j1ki1, ind6, T(ind6), N, N);
        A7 = sparse(jk1i1, ind7, T(ind7), N, N);
        A8 = sparse(j1k1i1, ind8, T(ind8), N, N);

        v1dxz = zeros(N, 1); v2dxz = zeros(N, 1);
        v3dxz = zeros(N, 1); v4dxz = zeros(N, 1);
        v5dxdz = zeros(N, 1); v6dxdz = zeros(N, 1);
        v7dxdz = zeros(N, 1); v8dxdz = zeros(N, 1);

        v1dyz = zeros(N, 1); v2dyz = zeros(N, 1);
        v3dyz = zeros(N, 1); v4dyz = zeros(N, 1);
        v5dydz = zeros(N, 1); v6dydz = zeros(N, 1);
        v7dydz = zeros(N, 1); v8dydz = zeros(N, 1);

        v1dzz = zeros(N, 1); v2dzz = zeros(N, 1);
        v3dzz = zeros(N, 1); v4dzz = zeros(N, 1);
        v5dzdz = zeros(N, 1); v6dzdz = zeros(N, 1);
        v7dzdz = zeros(N, 1); v8dzdz = zeros(N, 1);

        v1dxz(ind1) = (-1) * (1 - yp(ind1)) .* (1 - zp(ind1));
        v2dxz(ind2) = (1 - yp(ind2)) .* (1 - zp(ind2));
        v3dxz(ind3) = (-1) * (yp(ind3)) .* (1 - zp(ind3));
        v4dxz(ind4) = (yp(ind4)) .* (1 - zp(ind4));
        v5dxdz(ind5) = (-1) * (1 - yp(ind5)) .* (zp(ind5));
        v6dxdz(ind6) = (1 - yp(ind6)) .* (zp(ind6));
        v7dxdz(ind7) = (-1) * (yp(ind7)) .* (zp(ind7));
        v8dxdz(ind8) = (yp(ind8)) .* (zp(ind8));

        Txx = A1 * sdiag(v1dxz) + A2 * sdiag(v2dxz) + A3 * sdiag(v3dxz) + A4 * sdiag(v4dxz) + ...
            A5 * sdiag(v5dxdz) + A6 * sdiag(v6dxdz) + A7 * sdiag(v7dxdz) + A8 * sdiag(v8dxdz);

        v1dyz(ind1) = (1 - xp(ind1)) .* (-1) .* (1 - zp(ind1));
        v2dyz(ind2) = (xp(ind2)) .* (-1) .* (1 - zp(ind2));
        v3dyz(ind3) = (1 - xp(ind3)) .* (1 - zp(ind3));
        v4dyz(ind4) = xp(ind4) .* (1 - zp(ind4));
        v5dydz(ind5) = (1 - xp(ind5)) .* (-1) .* (zp(ind5));
        v6dydz(ind6) = (xp(ind6)) .* (-1) .* (zp(ind6));
        v7dydz(ind7) = (1 - xp(ind7)) .* (zp(ind7));
        v8dydz(ind8) = xp(ind8) .* zp(ind8);

        Tyy = A1 * sdiag(v1dyz) + A2 * sdiag(v2dyz) + A3 * sdiag(v3dyz) + A4 * sdiag(v4dyz) + ...
            A5 * sdiag(v5dydz) + A6 * sdiag(v6dydz) + A7 * sdiag(v7dydz) + A8 * sdiag(v8dydz);

        v1dzz(ind1) = (1 - xp(ind1)) .* (1 - yp(ind1)) .* (-1);
        v2dzz(ind2) = xp(ind2) .* (1 - yp(ind2)) .* (-1);
        v3dzz(ind3) = (1 - xp(ind3)) .* (yp(ind3)) .* (-1);
        v4dzz(ind4) = (xp(ind4)) .* (yp(ind4)) .* (-1);
        v5dzdz(ind5) = (1 - xp(ind5)) .* (1 - yp(ind5));
        v6dzdz(ind6) = xp(ind6) .* (1 - yp(ind6));
        v7dzdz(ind7) = (1 - xp(ind7)) .* yp(ind7);
        v8dzdz(ind8) = xp(ind8) .* yp(ind8);

        Tzz = A1 * sdiag(v1dzz) + A2 * sdiag(v2dzz) + A3 * sdiag(v3dzz) + A4 * sdiag(v4dzz) + ...
            A5 * sdiag(v5dzdz) + A6 * sdiag(v6dzdz) + A7 * sdiag(v7dzdz) + A8 * sdiag(v8dzdz);

        Tx = Txx / hx;
        Ty = Tyy / hy;
        Tz = Tzz / hz;
    end

end
