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

############
def prime():

#  '2-3-512' : {'lA' : 2, 'lB' : 3, 'eA' : 258, 'eB' : 161, 'f' : 186, 'pm1' : -1}
    la = 2
    lb = 3
    ea = 130
    eb = 81
    f = 22
    pm1 = -1
    
    rangeLim = 100
        
        
    for a in range(0,rangeLim):
      
        ea = ea-1
        p = "%s^%s * %s^%s" % (la, ea, lb, eb) + (f != 1)*(" * %s" % f) + (pm1 == 1)*" + 1" + (pm1 == -1)*" - 1"
        p = sage_eval(p)
        
        if is_prime(p):
            print "%s-bits prime p = %s" % (p.nbits(), p)
            print "ea: %d" % (ea)
            print "eb: %d" % (eb)




        eb = eb-1
        p = "%s^%s * %s^%s" % (la, ea, lb, eb) + (f != 1)*(" * %s" % f) + (pm1 == 1)*" + 1" + (pm1 == -1)*" - 1"
        p = sage_eval(p)
            
        if is_prime(p):
            print "%s-bits prime p = %s" % (p.nbits(), p)
            print "ea: %d" % (ea)
            print "eb: %d" % (eb)

        


    