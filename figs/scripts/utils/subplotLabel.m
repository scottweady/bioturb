function h = subplotLabel(s, ax, varargin)

    h = annotation('textbox', 'units', 'inches', 'fontsize', 18, 'string', s, 'interpreter', 'latex', 'edgecolor', 'none');

    if nargin > 2; pos = varargin{1};
    else; pos = 'northwest';
    end

    if nargin > 3; offset = varargin{2};
    else; offset = [0 0];
    end

    switch pos
        case 'northeast'
            h.Position(1) = ax.Position(1) + ax.Position(3);
            h.Position(2) = ax.Position(2) + ax.Position(4);
        case 'southeast'
            h.Position(1) = ax.Position(1) + ax.Position(3);
            h.Position(2) = ax.Position(2);
        case 'southwest'
            h.Position(1) = ax.Position(1);
            h.Position(2) = ax.Position(2);
        case 'northwest'
            h.Position(1) = ax.Position(1);
            h.Position(2) = ax.Position(2) + ax.Position(4);
    end

    h.Position(1:2) = h.Position(1:2) + offset;
