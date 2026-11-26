"""
Exact geometry of one step: crossing of forbidden cells, entry into the goal circle, map borders (sliding with friction).

Cells are the unit squares [i, i+1) x [j, j+1). A segment enters a cell when it goes through the inside of the cell;
touching a forbidden cell only at a corner does not count. A segment that ends exactly on the lower or left edge
of a forbidden cell does count, because its end point belongs to that cell.
"""
import math


def segment_first_blocked(p0, p1, blocked):
    """
    Fraction s in [0, 1] of the segment p0 -> p1 at which it first enters a cell where blocked[y, x] is True,
    or None. The cell containing p0 is not tested (the boat is never inside a forbidden cell).
    Grid traversal of Amanatides and Woo (1987): every cell crossed by the segment is visited, in order.
    """
    x0, y0 = float(p0[0]), float(p0[1])
    x1, y1 = float(p1[0]), float(p1[1])
    ix, iy = math.floor(x0), math.floor(y0)
    dx, dy = x1 - x0, y1 - y0
    inf = math.inf
    if dx > 0:
        sx, tmx, tdx = 1, (ix + 1 - x0) / dx, 1.0 / dx
    elif dx < 0:
        sx, tmx, tdx = -1, (ix - x0) / dx, -1.0 / dx
    else:
        sx, tmx, tdx = 0, inf, inf
    if dy > 0:
        sy, tmy, tdy = 1, (iy + 1 - y0) / dy, 1.0 / dy
    elif dy < 0:
        sy, tmy, tdy = -1, (iy - y0) / dy, -1.0 / dy
    else:
        sy, tmy, tdy = 0, inf, inf
    height, width = blocked.shape
    while True:
        t = tmx if tmx < tmy else tmy
        if t >= 1.0:
            break
        if tmx < tmy:
            ix += sx
            tmx += tdx
        elif tmy < tmx:
            iy += sy
            tmy += tdy
        else:  # the segment goes exactly through a corner: it does not enter the two side cells
            ix += sx
            iy += sy
            tmx += tdx
            tmy += tdy
        if not (0 <= ix < width and 0 <= iy < height):
            break
        if blocked[iy, ix]:
            return t if t > 0.0 else 0.0
    # the end point may lie exactly on the edge of the next cell: it belongs to that cell
    ex, ey = math.floor(x1), math.floor(y1)
    if (ex, ey) != (ix, iy) and 0 <= ex < width and 0 <= ey < height and blocked[ey, ex]:
        return 1.0
    return None


def segment_circle_entry(p0, p1, center, radius):
    """Fraction s in [0, 1] of the segment p0 -> p1 at which it enters the closed disc (center, radius), or None."""
    fx, fy = float(p0[0]) - center[0], float(p0[1]) - center[1]
    dx, dy = float(p1[0]) - float(p0[0]), float(p1[1]) - float(p0[1])
    c = fx * fx + fy * fy - radius * radius
    if c <= 0.0:
        return 0.0  # already inside
    a = dx * dx + dy * dy
    b = 2.0 * (fx * dx + fy * dy)
    if a == 0.0 or b >= 0.0:  # not moving, or not moving towards the centre: no entry
        return None
    disc = b * b - 4.0 * a * c
    if disc < 0.0:
        return None
    s = 2.0 * c / (-b + math.sqrt(disc))  # smaller root, numerically stable form
    return s if s <= 1.0 else None


def apply_borders(x0, y0, vx, vy, size, friction):
    """
    Move from (x0, y0) by (vx, vy), staying inside the map [0, size)^2. Returns (x, y, vx, vy).
    A step that would leave the map across one border stops on that border (that velocity component becomes zero)
    and slides along it with friction: the other component, of the velocity and of this step's move, is multiplied
    by `friction`. A step that would leave across two borders at once stops in the corner.
    """
    top = math.nextafter(size, 0.0)
    px, py = x0 + vx, y0 + vy
    out_x = not (0.0 <= px < size)
    out_y = not (0.0 <= py < size)
    if out_x:
        px = 0.0 if px < 0.0 else top
    if out_y:
        py = 0.0 if py < 0.0 else top
    if out_x and out_y:
        return px, py, 0.0, 0.0
    if out_x:
        vy *= friction
        return px, y0 + vy, 0.0, vy
    if out_y:
        vx *= friction
        return x0 + vx, py, vx, 0.0
    return px, py, vx, vy
