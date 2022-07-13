# Hydrophobic Moment Calculator

Script to calculate a set of properties from a protein sequence:
  - hydrophobicity (according to a particular scale)
  - mean hydrophobic dipole moment assuming it is an alpha-helix.
  - total charge (at pH 7.4)
  - amino acid composition
  - discimination factor according to Rob Keller (IJMS, 2011)

Reproduces the same values as the HeliQuest server, but it's a script.

Copied from [my own gist](https://gist.github.com/JoaoRodrigues/568c845915aea3efa3578babfd72423c)
to make it simpler for others to re-use.

## Contact
Joao Rodrigues (j.p.g.l.m.rodrigues@gmail.com)

## Example Usage
```
# Read sequences from a file
python hydrophobic_moment.py -f example_data/example.fasta

# Read single sequence from command-line
python hydrophobic_moment.py -s AAASSDADASDASEAWD
```
