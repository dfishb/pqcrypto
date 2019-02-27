# Copyright (c) 2011 Luca De Feo.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from sage.misc.sage_timeit import sage_timeit
import sage.misc.misc

la=2
lb=3

ea=121
eb=1
f=1
t=-1
k=0
lower = 2^224
upper = 2^256


def param_search():

    for ea in range(110,150):
        print("\nea: %d",ea);

        for eb in range(1,230):

            for f in range(1,230):
        
                num = (la^ea)*(lb^eb)*f+1
                

                if is_prime(num):
                    k=0
                    while 2**k < num:
                        k = k+1
                     #   print k
                        
                        if k==240 & num>lower & num<upper:
                            print("la: %d\nlb: %d\nea: %d\neb: %d\nlf: %d\nt: %d\n",la,lb,ea,eb,f,t);
                  #  ss_isogeny_gen(lA, lB, eA, eB, f, pm1)
                       
                         
                num = (la^ea)*(lb^eb)*f-1
                if is_prime(num):
                       k=0
                       while 2**k < num: 
                            k = k+1
                         #   print k
                            
                       if k==240 & num>lower & num<upper :
                            print("la: %d\nlb: %d\nea: %d\neb: %d\nlf: %d\nt: %d\n",la,lb,ea,eb,f,t);
                  #  ss_isogeny_gen(lA, lB, eA, eB, f, pm1)
                  
                  
                  