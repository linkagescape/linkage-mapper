##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: gapdt.py 322 2008-11-29 19:57:22Z viral $
##

from numpy import *
from scipy import sparse

class gapdt:

    def sprand (self, x):
        s = sparse.coo_matrix(x);

    def sprand (self, m, n, density):        
        nnz = fix(m*n*density)
        i = random.random_integers(0, m-1, nnz)
        j = random.random_integers(0, n-1, nnz)
        ij = c_[i, j].T
        data = random.rand(nnz)
        
        s = sparse.coo_matrix((data, ij), shape=(m, n))
        return s

    # def relabel(self, oldlabel, offset=0):
        # newlabel = zeros(size(oldlabel), dtype='int32')
        # s = sort(oldlabel)
        # perm = argsort(oldlabel)
        # f = where(diff(concatenate(([s[0]-1], s))))
        # newlabel[f] = 1
        # newlabel = cumsum(newlabel)
        # newlabel[perm] = copy(newlabel)
        # return newlabel-1+offset

                
    # def conditional_hooking (self, D, u, v):
        # Du = D[u]
        # Dv = D[v]
        
        # hook = where ((Du == D[Du]) & (Dv < Du))
        # Du = Du[hook]
        # Dv = Dv[hook]
        
        # D[Du] = Dv
        # return D

    # def check_stars (self, D):
        # n = D.size
        # star = ones (n, dtype='int32')

        # notstars = where (D != D[D])
        # star[notstars] = 0
        # Dnotstars = D[notstars]
        # star[Dnotstars] = 0
        # star[D[Dnotstars]] = 0

        # star = star[D]
        # return star

    # def pointer_jumping (self, D):
        # n = D.size
        # Dold = zeros(n, dtype='int32');

        # while any(Dold != D):
            # Dold = D
            # D = D[D]

        # return D

    # def components(self, G):
        # G = sparse.coo_matrix(G)
        # n = G.shape[0]
        # D = arange (0, n, dtype='int32')
        # U = G.row
        # V = G.col

        # while True:
            # D = self.conditional_hooking(D, U, V)
            # star = self.check_stars (D)
            
            # if (sum(star) == n):
                # return self.relabel(D, 1)
                # break

            # D = self.pointer_jumping(D)

    def relabel(self, oldlabel, offset=0):
        newlabel = zeros(size(oldlabel), dtype='int32')
        s = sort(oldlabel)
        perm = argsort(oldlabel)
        f = where(diff(concatenate(([s[0]-1], s))))
        newlabel[f] = 1
        newlabel = cumsum(newlabel)
        newlabel[perm] = copy(newlabel)
        return newlabel-1+offset

                
    def conditional_hooking (self, D, star, u, v):
        Du = D[u]
        Dv = D[v]
        
        hook = where ((star[u] == 1) & (Du > Dv))
        D[Du[hook]] = Dv[hook]

        return D

    def unconditional_hooking (self, D, star, u, v):
        Du = D[u]
        Dv = D[v]
        
        hook = where((star[u] == 1) & (Du != Dv))
        D[Du[hook]] = Dv[hook]

        return D

    def check_stars (self, D, star):
        star[:] = 1
        notstars = where (D != D[D])
        star[notstars] = 0
        star[D[D[notstars]]] = 0
        star = star[D]
        return star

    def components(self, G):
        Gcoo = sparse.coo_matrix(G)
        n = G.shape[0]
        U = Gcoo.row
        V = Gcoo.col

        D = arange (0, n, dtype='int32')
        star = zeros(n, 'int32')

        all_stars = False
        while True:
            star = self.check_stars (D, star)
            D = self.conditional_hooking(D, star, U, V)

            star = self.check_stars (D, star)
            D = self.unconditional_hooking (D, star, U, V)

            # pointer jumping
            D = D[D]

            if all_stars == True:
                return self.relabel(D, 1)
            
            if sum(star) == n:
                all_stars = True

        
    def subsref(self, A, I, J):
        B = A[:, J][I, :]
        
        return B

    def deleterowcol(self, A, delrow, delcol):
        m = A.shape[0]
        n = A.shape[1]

        keeprows = delete (arange(0, m), delrow)
        keepcols = delete (arange(0, n), delcol)

        return A[keeprows][:,keepcols]
