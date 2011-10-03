      subroutine simark(n, m, sumQ, runif, curr, x)
      implicit none
      integer n, m, j, k, curr
      integer x(n)
      double precision runif(n)
      double precision sumQ(m,m)
c     G. Filion: n is the length of the sequence
c     G. Filion: m is the dimension of the matrix
c     G. Filion: sumQ contains cumulative transition probabilities
c     G. Filion: runif is the vector of simulated probabilities
c     G. Filion: curr is the current state
c     G. Filion: x is the sequence
      j = 2
      do while(j .le. n)
           k = 1
           do while(k .le. m)
               if ( runif(j) .le. sumQ(curr,k) ) then
                   curr = k
                   x(j) = k
                   k = m
               endif
               k = k+1
           enddo
           j = j+1
      enddo
      end
