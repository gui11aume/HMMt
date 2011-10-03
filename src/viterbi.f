      subroutine vit(n, m, logini, logPem, logQ, path, maxmat)
c     written by Guillaume Filion, 05.24.2008
      implicit none
      integer i, j, k, n, m
      double precision tmp, themax, logini(m)
      double precision logPem(n,m), logQ(m,m), maxmat(n,m)
      integer path(n)
c     G. Filion: the first forward step is handled separately
      j = 1
      do while(j .le. m)
          maxmat(1,j) = logini(j) + logPem(1,j)
          j = j+1
      enddo
      i = 2
      do while(i .le. n)
          j = 1
          do while(j .le. m)
              themax = maxmat(i-1,1) + logQ(1,j)
              k = 2
              do while(k .le. m)
                  tmp = maxmat(i-1,k) + logQ(k,j)
                  if (tmp .gt. themax) themax = tmp
                  k = k+1
              enddo
              maxmat(i,j) = themax + logPem(i,j)
              j = j+1
          enddo
          i = i+1
      enddo
c     G. Filion: start of the backward loop
c     G. Filion: the first step is handled separately
      k = 1
      j = 2
      do while(j .le. m)
          if (maxmat(n,j) .gt. maxmat(n,k)) k = j
          j = j+1
      enddo
      path(n) = k
      i = n-1
      do while(i .ge. 1)
      themax = logQ(1, path(i+1)) + maxmat(i,1)
      k = 1
      j = 2
          do while(j .le. m)
              tmp = logQ(j, path(i+1)) + maxmat(i,j)
              if (tmp .gt. themax) then
                  themax = tmp
                  k = j
              endif
              j = j+1
          enddo
          path(i) = k
          i = i-1
      enddo
      end
