#!/usr/bin/env python
"""Share the products that 81 observations write once per system.

    env/bin/python dedup_image_products.py [--apply]

An observation that covers more than one system of the sample is imaged and
detected once per system, into each system's own directory, so its products
exist several times over. They are the same products: the detection runs on the
whole field, not on a sub-region around the system, so two copies of a source
list differ only in the DATE and CHECKSUM cards recording when each was
written. Checked on obsid 144, whose two copies have identical data in every
HDU and identical `source_list.csv` and `src.reg`.

Nothing is deleted here. Copies whose *data* are identical are replaced by hard
links to one inode, so every path still exists and every system directory still
lists everything it did before; only the bytes behind them are shared. A group
whose data differ is left alone and reported.

Without --apply this only measures.
"""

import os
import re
import sys
from collections import defaultdict

import numpy as np
from astropy.io import fits

base = "/home/marina/Doctorado/2026/HMXB_project/Observations/Images"
APPLY = "--apply" in sys.argv
GB = 1024.0 ** 3


def obsid_of(subdir, name):
    if subdir == "Chandra_images":
        m = re.match(r"^(\d+)_", name)
    else:
        m = re.search(r"(?:^|_)(0\d{9})[_.]", name)
    return m.group(1) if m else None


def equal(x, y):
    """np.array_equal, but NaN equals NaN: a detection table is full of them
    and NaN != NaN would report every float column as different."""
    try:
        return np.array_equal(x, y, equal_nan=True)
    except TypeError:
        return np.array_equal(x, y)


def same_data(a, b):
    """Equal in every data unit, whatever the headers say."""
    if os.path.getsize(a) != os.path.getsize(b):
        return False
    if not a.endswith((".fits", ".fits.gz")):
        return open(a, "rb").read() == open(b, "rb").read()
    try:
        with fits.open(a, memmap=True) as ha, fits.open(b, memmap=True) as hb:
            if len(ha) != len(hb):
                return False
            for x, y in zip(ha, hb):
                if x.data is None and y.data is None:
                    continue
                if x.data is None or y.data is None:
                    return False
                if getattr(x.data, "dtype", None) is not None and x.data.dtype.names:
                    if x.data.dtype.names != y.data.dtype.names:
                        return False
                    if not all(equal(x.data[n], y.data[n])
                               for n in x.data.dtype.names):
                        return False
                elif not equal(x.data, y.data):
                    return False
        return True
    except Exception:
        return False


groups = defaultdict(list)
for system in sorted(os.listdir(base)):
    for subdir in ("Chandra_images", "XMM-Newton_images"):
        d = os.path.join(base, system, subdir)
        if not os.path.isdir(d):
            continue
        for name in os.listdir(d):
            o = obsid_of(subdir, name)
            if o:
                groups[(subdir, o, name)].append(os.path.join(d, name))

shared = reclaimed = 0
already = 0
differ = []
for key, paths in sorted(groups.items()):
    if len(paths) < 2:
        continue
    keep, rest = paths[0], paths[1:]
    for p in rest:
        if os.stat(p).st_ino == os.stat(keep).st_ino:
            already += 1
            continue
        if not same_data(keep, p):
            differ.append((key[1], key[2], keep, p))
            continue
        size = os.path.getsize(p)
        if APPLY:
            tmp = p + ".relink"
            os.link(keep, tmp)
            os.replace(tmp, p)
        shared += 1
        reclaimed += size

print("%s" % ("applied" if APPLY else "dry run, nothing changed"))
print("  files already sharing an inode : %d" % already)
print("  files linked to an existing copy: %d" % shared)
print("  bytes no longer duplicated      : %.1f GB" % (reclaimed / GB))
if differ:
    print("  groups left alone, data differ  : %d" % len(differ))
    for o, name, a, b in differ[:20]:
        print("    %-10s %-34s %s | %s" % (o, name, a.split("/")[-3], b.split("/")[-3]))
