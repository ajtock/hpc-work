#!/usr/bin/env python

from time import time

start = time()
print("this")
end = time()
print("Elapsed: %.19f" % (end - start) + "s")
print("Elapsed: %s" % (end - start) + "s")
print(f"Elapsed: {end - start:.19f}s")
