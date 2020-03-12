c Copyright 2005-2020 Kendall F. Morris

c This file is part of the Xanalysis software suite.

c     The Xanalysis software suite is free software: you can redistribute
c     it and/or modify it under the terms of the GNU General Public
c     License as published by the Free Software Foundation, either
c     version 3 of the License, or (at your option) any later version.

c     The suite is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.

c     You should have received a copy of the GNU General Public License
c     along with the suite.  If not, see <https://www.gnu.org/licenses/>.

      module mod_mean_and_sd_routines
      contains
c       filename: mean_and_sd_routines          lss
c       created: 04-Dec-2003            lss

c       These subroutines will calculate the mean and standard deviation of data contained within a one-dimensional array of
c       double-precision floats, single-precision floats, or integers.

c       The calling code must supply the length of the array as the second argument.


      subroutine dp_real_mean_and_sd(data_arr,length,mean,sd)

      double precision data_arr(length),mean,sd,sum_sqr,sum
      integer length

      mean = 0.0
      sum = 0.0
      do i = 1, length
         mean = mean + (data_arr(i)/float(length))
         sum = sum + data_arr(i)
      end do
      mean = sum/float(length)

      sum_sqr = 0.0
      do i = 1, length
         sum_sqr = sum_sqr + ((mean - data_arr(i))**2)
      end do
      sd = SQRT(sum_sqr/(float(length)-1.))
      
      return
      end

*     <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
*     <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

      subroutine sp_real_mean_and_sd(data_arr,length,mean,sd)

      real data_arr(length),mean,sd,sum_sqr,sum
      integer length

      mean = 0.0
      sum = 0.0
      do i = 1, length
         mean = mean + (data_arr(i)/float(length))
         sum = sum + data_arr(i)
      end do
      mean = sum/float(length)

      sum_sqr = 0.0
      do i = 1, length
         sum_sqr = sum_sqr + ((mean - data_arr(i))**2)
      end do
      sd = SQRT(sum_sqr/(float(length)-1.))
      
      return
      end
      
*     <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
*     <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

      subroutine int_mean_and_sd(data_arr,length,mean,sd)

      real sd, sum_sqr
      integer*4 length, mean, data_arr(length), sum

      mean = 0
      sum = 0
      do i = 1, length
         sum = sum + data_arr(i)
      end do
      mean = NINT(float(sum)/float(length))
      mean = sum/length

      sum_sqr = 0.0
      do i = 1, length
         sum_sqr = sum_sqr + float((mean - data_arr(i))**2)
      end do
      sd = SQRT(sum_sqr/(float(length)-1.))
      
      return
      end
      
      end module mod_mean_and_sd_routines
