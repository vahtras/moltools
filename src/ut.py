def upper_triangular(n, start=0):
    """Recursive generator for triangular looping of Carteesian tensor"""
    if n > 2:
        for i in range(start, 3):
            for j in upper_triangular(n-1, start=i):
                yield (i,) + j
    else:
        for i in range(start, 3):
            for j in range(i, 3):
                yield i, j
