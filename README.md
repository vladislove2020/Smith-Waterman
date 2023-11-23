##Smith-Waterman algorithm implementation

Smith-Waterman algorithm is a popular DynProg algorithm to find the optimal local alignment of two strings.
It is often used in bioinformatics to discover local similarities in protein or nucleic acid sequences.
It can also be used to work with languages.
More info about the algorithm in Wikipedia https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

NOTE 1:
This implementation is educational and very simple and doesn't take into account probabilities of transitions between symbols.
So it cannot be useful to imitate true evolutionary alignment of protein sequences, because, for example, transition Asp -> Trp is many times less probable than transition Asp -> Glu.
If you wanna try professional Smith-Waterman algorithm implementation for your bioinformatics tasks, you can use, for example, "water" program from EMBOSS package.
If you wanna just have fun and see how SW algorithm works on easy examples, this repository can be useful :)

NOTE2:
The program implements two functions for aligning - with linear and affine gap penalties.
Penalty is a decrease of alignment score (SW algorithm searches for the alignment with maximum score).
Gap is a situation, when some symbols in sequence 1 don't match any symbols in sequence 2.
In bioinformatics it usually means that unmatched symbols were either removed from sequence 2 or inserted into sequence 1 during the evolutionary process.
If you use linear gap penalties, 5 one-letter gaps and 1 five-letters gap give you the same penalty.
If you use affine gap penalties, opening the gap is penalted more than prolongating the gap. 
Affine gap penalties provide more correct approach, when you align biological sequences, because few long insertions/deletions are more probable than many short ones.
