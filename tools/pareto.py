"""pareto: non-dominated filtering and shortlist selection for two minimised objectives."""
import numpy as np


def pareto_indices(x, y) -> list[int]:
    """Indices of non-dominated points when minimising both x and y, sorted by x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    order = np.lexsort((y, x))
    front, best_y = [], np.inf
    for i in order:
        if y[i] < best_y:
            front.append(int(i))
            best_y = y[i]
    return front


def spread_along_front(x_front, n: int) -> list[int]:
    """Pick up to n positions along a front (sorted by x), evenly spaced in log10(x), endpoints included."""
    logx = np.log10(np.asarray(x_front, float))
    if len(logx) <= n:
        return list(range(len(logx)))
    targets = np.linspace(logx[0], logx[-1], n)
    picks = []
    for target in targets:
        for j in np.argsort(np.abs(logx - target)):
            if int(j) not in picks:
                picks.append(int(j))
                break
    return sorted(picks)


def knee_index(x_front, y_front) -> int:
    """Knee of a front: the point farthest from the chord between its endpoints, in log-log space."""
    lx, ly = np.log10(np.asarray(x_front, float)), np.log10(np.asarray(y_front, float))
    if len(lx) < 3:
        return 0
    dx, dy = lx[-1] - lx[0], ly[-1] - ly[0]
    dist = np.abs(dx * (ly - ly[0]) - dy * (lx - lx[0])) / np.hypot(dx, dy)
    return int(np.argmax(dist))
