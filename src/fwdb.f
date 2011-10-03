      subroutine fwdb(m, n, cur, prob, Q, ck, phi, trans, sumphi)
c     modified from the package HiddenMarkov by Guillaume Filion, 04.27.2008
c     G. Filion: notations are borrowed from Inference in Hidden Markov
      implicit none
      integer i, j, k, m, n
      double precision sumphi, cur(m)
      double precision prob(n,m), Q(m,m)
      double precision ck(n), phi(n,m)
      double precision back(m,m)
      double precision trans(m,m)
      double precision tmp(m)
c     G. Filion: cur is the working vector (initialized with initialProb)
c     G. Filion: the algorithm replaces prob by the forward alphas
c     G. Filion: the algorithm does not compute the betas
c     G. Filion: instead, it updates the transition probabilities
c     G. Filion: the first forward step is handled separately
      j = 1
      do while(j .le. m)
           cur(j) = cur(j)*prob(1,j)
           ck(1) = ck(1) + cur(j)
           j = j+1
      enddo
      j = 1
      do while(j .le. m)
           cur(j) = cur(j)/ck(1)
c     G. Filion: replace prob by the alphas
           prob(1,j) = cur(j)
           j = j+1
      enddo
c     G. Filion: start of the forward loop
      i = 2
      do while(i .le. n)
c     G. Filion: forward transitions
          j = 1
          do while(j .le. m)
              k = 1
              tmp(j) = 0.0
              do while(k .le. m)
                  tmp(j) = tmp(j) + cur(k)*Q(k,j)
                  k = k+1
              enddo
              j = j+1
          enddo
          j = 1
          do while(j .le. m)
              cur(j) = tmp(j)
              j = j+1
          enddo
c     G. Filion: normalization
          j = 1
          do while(j .le. m)
              cur(j) = cur(j)*prob(i,j)
              ck(i) = ck(i) + cur(j)
              j = j+1
          enddo
          j = 1
          do while(j .le. m)
              cur(j) = cur(j)/ck(i)
c     G. Filion: replacing prob by the alphas
              prob(i,j) = cur(j)
              j = j+1
          enddo
          i = i+1
      enddo
c     G. Filion: the first backward step is handled separately
c     G. Filion: prob now contains the alpha forward probabilities
      j = 1
      do while(j .le. m)
          phi(n,j) = prob(n,j)
          j = j+1
      enddo
c     G. Filion: start of the backward loop
      i = n-1
      do while(i .ge. 1)
c     G. Filion: backward transitions (and normaliation)
          j = 1
          do while(j .le. m)
              k = 1
              sumphi = 0.0
              do while(k .le. m)
                 back(j,k) = prob(i,k)*Q(k,j)
                 sumphi = sumphi + back(j,k)
                 k = k+1
              enddo
              k = 1
              do while(k .le. m)
                 back(j,k) = back(j,k) / sumphi
                 k = k+1
              enddo
              j = j+1
          enddo
          j = 1
c     G. Filion: computing the posterior probabilites and transitions
          do while(j .le. m)
              k = 1
              do while(k .le. m)
c     G. Filion: directly summing the posterior transitions
                  trans(j,k) = trans(j,k) + phi(i+1,k)*back(k,j)
                  phi(i,j) = phi(i,j) + phi(i+1,k)*back(k,j)
                  k = k+1
              enddo
          j = j+1
          enddo
          i = i-1
      enddo
c     G. Filion: I commented out this step because it is now carried out in R
c     G. Filion: normalizing the posterior transition probabilities
c     G. Filion: NB: sumphi is used as a temporary variable
c      j = 1
c      do while(j .le. m)
c          k = 1
c          sumphi = 0.0
c          do while(k .le. m)
c               sumphi = sumphi + trans(j,k)
c               k = k+1
c          enddo
c          k = 1
c          do while(k .le. m)
c               trans(j,k) = trans(j,k) / sumphi
c               k = k+1
c          enddo
c          j = j+1
c      enddo
c     G. Filion: storing the log-likelihodd in sumphi
      sumphi = log(ck(1))
      i = 2
      do while(i .le. n)
          sumphi = sumphi + log(ck(i))
          i = i+1
      enddo
      end
