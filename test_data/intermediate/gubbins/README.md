# TEST gubbins output for stage 8

NOT BIOLOGICAL RESULTS and not gubbins output from any real run.
Committed so TEST can exercise stage 8's real parser and writer with
no external tool installed, exactly as `stages.similarity.run` reads a
fixture.

Regenerate with
`python3 artifacts/round11/dag/make_test_gubbins_fixture.py`
(deterministic - same bytes every time).

Files, under the names `adapters.gubbins` looks for:

| file | suffix constant |
| --- | --- |
| `recombination.per_branch_statistics.csv` | `PER_BRANCH_STATISTICS_SUFFIX` |
| `recombination.node_labelled.final_tree.tre` | `NODE_LABELLED_TREE_SUFFIX` |
| `recombination.filtered_polymorphic_sites.fasta` | `SNP_ALIGNMENT_SUFFIX` |

The tree labels internal nodes after the closing parenthesis, which is
gubbins' shape and not Newick. Two tips carry zero SNPs on purpose.
