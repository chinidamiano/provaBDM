import pandas as pd

their_sigs = pd.read_pickle("./paper_signature.p")
our_sigs = pd.read_pickle("./output/1-Replication/signature.p")

their_sigs = their_sigs.to_dict()

tmp = []
for key in their_sigs["coef"]:
    tmp.append(key[1])

intersection = list(set.intersection(set(tmp), set(our_sigs.index)))
intersection.sort()

with open("./output/1-Replication/signature_intersection.txt", 'w') as handle:
    handle.write("\n".join(intersection))
