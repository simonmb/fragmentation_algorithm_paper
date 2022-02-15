## fragmentation_algorithm_paper
this is the immutable version of the algorithm for the following paper:
- [**Flexible Heuristic Algorithm for Automatic Molecule Fragmentation: Application to the UNIFAC Group Contribution Model**](https://doi.org/10.1186/s13321-019-0382-3)

# ⚠️If you are interested in a new version⚠️
- [New version of the fragmentation algorithm](https://github.com/simonmb/fragmentation_algorithm)
- Why is there a newer version?
  - The fragmentation algorithm originally published with [the paper](https://doi.org/10.1186/s13321-019-0382-3) tried to find smaller functional groups that are contained within other larger functional groups to sort them automatically in an intelligent way. This turned out to be quite difficult to implement as the capabilities of RDKit to match SMARTS with SMARTS are  limited. This lead me to writing workarounds that became broken in subsequent RDKit versions.
  - The newer version does not try to automatically sort the groups you are searching for but relies on the user to provide this order. This allows the algorithm to be applicable with recent RDKit versions without problems.
