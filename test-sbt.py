import sourmash

## build an sbt
sbt = sourmash.create_sbt_index()
mh = sourmash.MinHash(ksize=7, n=0, scaled=1, track_abundance=True, is_protein=True)
mh.add_sequence("MKVIRVG")
sig = sourmash.SourmashSignature(mh, name="test")
print("md5: " + sig.md5sum())
print(mh.get_mins())
leaf = sourmash.sbtmh.SigLeaf(sig.md5sum(), sig)
sbt.add_node(leaf)
sbt.save("test.sbt.zip")

# everything the same
sbt = sourmash.create_sbt_index()
mh = sourmash.MinHash(ksize=7, n=0, scaled=1, track_abundance=True, is_protein=True)
mh.add_sequence("MKVIRVG")
sig = sourmash.SourmashSignature(mh, name="test")
print("md5: " + sig.md5sum())
print(sig.minhash.get_mins())
leaf = sourmash.sbtmh.SigLeaf(sig.md5sum(), sig)
sbt.add_node(leaf)
sbt.save("test.sbt.zip")

## different sequence (same name)
sbt = sourmash.create_sbt_index()
mh = sourmash.MinHash(ksize=7, n=0, scaled=1, track_abundance=True, is_protein=True)
mh.add_sequence("VGLGDRG")
sig = sourmash.SourmashSignature(mh, name="test")
print("md5: " + sig.md5sum())
print(mh.get_mins())
leaf = sourmash.sbtmh.SigLeaf(sig.md5sum(), sig)
sbt.add_node(leaf)
sbt.save("test.sbt.zip")

## compared to previous: same sequence sequence, different name
sbt = sourmash.create_sbt_index()
mh = sourmash.MinHash(ksize=7, n=0, scaled=1, track_abundance=True, is_protein=True)
mh.add_sequence("VGLGDRG")
sig = sourmash.SourmashSignature(mh, name="test2")
print("md5: " + sig.md5sum())
print(mh.get_mins())
leaf = sourmash.sbtmh.SigLeaf(sig.md5sum(), sig)
sbt.add_node(leaf)
sbt.save("test.sbt.zip")

# longer seq
sbt = sourmash.create_sbt_index()
mh = sourmash.MinHash(ksize=7, n=0, scaled=1, track_abundance=True, is_protein=True)
mh.add_sequence("MKVIRVGVVGLGDRGLH")
sig = sourmash.SourmashSignature(mh, name="test2")
print("md5: " + sig.md5sum())
print(sig.minhash.get_mins())
leaf = sourmash.sbtmh.SigLeaf(sig.md5sum(), sig)
sbt.add_node(leaf)
sbt.save("test.sbt.zip")

## compared to previous: same sequence and name (different ksize)
sbt = sourmash.create_sbt_index()
mh = sourmash.MinHash(ksize=6, n=0, scaled=1, track_abundance=True, is_protein=True)
mh.add_sequence("MKVIRVGVVGLGDRGLH")
sig = sourmash.SourmashSignature(mh, name="test2")
print("md5: " + sig.md5sum())
print(mh.get_mins())
leaf = sourmash.sbtmh.SigLeaf(sig.md5sum(), sig)
sbt.add_node(leaf)
sbt.save("test.sbt.zip")


