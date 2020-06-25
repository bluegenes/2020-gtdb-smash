#! /usr/bin/env python

# (from @ctb)

import sourmash
import sys
from sourmash.lca.lca_db import load_single_database

lcafile = sys.argv[1]
sbtfile = sys.argv[2]

db, _, _ = load_single_database(lcafile)

tree = sourmash.create_sbt_index()

for sig in db.signatures():
    tree.insert(sig)

tree.save(sbtfile)
