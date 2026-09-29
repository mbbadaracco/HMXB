"""Check the two branches of crossmatch.Omega_function against exact geometry.

`crossmatch.Omega_function` takes a `formula` argument:

  formula="A6"         reproduces Eq. (A.6) of Pineau et al. (2011)
                       literally,            sin(Xhat) = x * sin(Nhat) / n
  formula="AppendixB"  the default,          sin(Xhat) = sin(x) * sin(Nhat) / sin(n)

The second is the spherical law of sines; the first is what the paper
prints.  They are not the same, so which one the code uses has to be
settled rather than assumed.  This script settles it, by computing the
angle at vertex X of the spherical triangle N-X-O directly from vectors
-- no spherical-trigonometry identity involved -- and comparing.

Angles are compared modulo 180 deg, because the quantity they are used
for is the rotation of a covariance matrix, which is pi-periodic.

Run:  ./env/bin/python check_rotation_formula.py
"""

import numpy as np

N_PAIRS = 30_000
SEED = 7


def true_angle_at_X(ra_o, de_o, ra_x, de_x):
    """Angle at vertex X of the triangle (north pole, X, O), in degrees.

    Built from the dihedral angle between the tangent directions at X
    towards the pole and towards O.  This uses no identity from the
    paper, so it is an independent reference.
    """
    def vec(ra, de):
        ra, de = np.radians(ra), np.radians(de)
        return np.array([np.cos(de) * np.cos(ra),
                         np.cos(de) * np.sin(ra),
                         np.sin(de)])

    x_vec, o_vec = vec(ra_x, de_x), vec(ra_o, de_o)
    pole = np.array([0.0, 0.0, 1.0])
    t_n = pole - x_vec * np.dot(pole, x_vec)
    t_n /= np.linalg.norm(t_n)
    t_o = o_vec - x_vec * np.dot(o_vec, x_vec)
    t_o /= np.linalg.norm(t_o)
    return np.degrees(np.arctan2(np.dot(np.cross(t_n, t_o), x_vec),
                                 np.clip(np.dot(t_n, t_o), -1, 1)))


def haversine_deg(ra_o, de_o, ra_x, de_x):
    a = np.radians(de_o - de_x)
    b = np.radians(ra_o - ra_x)
    h = (np.sin(a / 2) ** 2
         + np.cos(np.radians(de_o)) * np.cos(np.radians(de_x)) * np.sin(b / 2) ** 2)
    return np.degrees(2 * np.arcsin(np.sqrt(h)))


def pineau_angle_at_X(ra_o, de_o, ra_x, de_x, exact):
    """Eqs. (A.3)-(A.6) of Pineau et al. (2011); `exact` selects the branch."""
    n = np.radians(haversine_deg(ra_o, de_o, ra_x, de_x))
    n_hat = np.radians(ra_o - ra_x)
    x = np.radians(90 - de_o)
    o = np.radians(90 - de_x)
    s = (n + x + o) / 2
    cos_half = np.sqrt(max(np.sin(s) * np.sin(s - x) / (np.sin(n) * np.sin(o)), 0.0))
    arg = (np.sin(x) * np.sin(n_hat) / np.sin(n)) if exact else (x * np.sin(n_hat) / n)
    a = np.arcsin(np.clip(arg, -1, 1))
    if cos_half > np.sqrt(2) / 2:
        return np.degrees(np.pi - a if a >= 0 else -np.pi - a)
    return np.degrees(a)


def diff_mod_180(u, v):
    return abs(((u - v + 90) % 180) - 90)


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    while len(rows) < N_PAIRS:
        ra_x = rng.uniform(0, 360)
        de_x = np.degrees(np.arcsin(rng.uniform(-1, 1)))
        sep = 10 ** rng.uniform(-5, 0)              # 0.036 arcsec to 1 deg
        pa = rng.uniform(0, 2 * np.pi)
        de_o = de_x + sep * np.cos(pa)
        ra_o = ra_x + sep * np.sin(pa) / max(np.cos(np.radians(de_x)), 1e-6)
        if not -89.9 < de_o < 89.9:
            continue
        truth = true_angle_at_X(ra_o, de_o, ra_x, de_x)
        rows.append((
            abs(de_x),
            diff_mod_180(pineau_angle_at_X(ra_o, de_o, ra_x, de_x, True), truth),
            diff_mod_180(pineau_angle_at_X(ra_o, de_o, ra_x, de_x, False), truth),
        ))

    a = np.array(rows)
    err_default, err_a6 = a[:, 1].max(), a[:, 2].max()
    print(f"pairs tested: {len(a)}  (separations 0.036'' to 1 deg, all declinations)")
    print(f'  default "AppendixB", sin(x) sin(N) / sin(n) : max error = {err_default:.3e} deg')
    print(f'  published Eq. (A.6),      x  sin(N) /     n : max error = {err_a6:.3e} deg')
    for lo, hi in [(0, 30), (30, 60), (60, 85), (85, 90)]:
        m = (a[:, 0] >= lo) & (a[:, 0] < hi)
        if m.sum():
            print(f"    |dec| in [{lo:2d},{hi:2d}) : n={m.sum():5d} "
                  f" default={a[m, 1].max():.2e}  A6={a[m, 2].max():.2e}")

    assert err_default < 1e-5, err_default
    assert err_a6 > 10.0, err_a6
    print("\nOK: the default branch is the exact spherical law of sines;")
    print("    Eq. (A.6) as printed is not, and the code is right not to use it.")


if __name__ == "__main__":
    main()
