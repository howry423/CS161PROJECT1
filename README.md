Run the following:

python3 python/solve_all.py inputs outputs1

This runs the sequential greedy algorithm on all inputs and stores outputs in the folder outputs1.

Then, modify the function solver() in solve_all.py as follows:

def solver(size: Size, instance: Instance) -> Solution:
    # Modify this function to use your imported solvers.
    # YOUR CODE HERE
    if size == Size.SMALL:
        return solve_naive(instance)
    elif size == Size.MEDIUM:
        return solve_naive(instance)
    elif size == Size.LARGE:
        return solve_naive(instance)


Then, run the following:

python3 python/solve_all.py inputs outputs2

This runs the naive algorithm on all inputs and stores outputs in the folder outputs2.

Lastly, run the following:

python3 python/merge.py --inputs inputs outputs1 outputs2 best

This merges the two output folders together into a folder in the root directory named "best". 