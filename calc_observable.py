#!/usr/local/bin/python3
import readline
import numpy as np
#
import Orbit
import Operator
import TransitionDensity

def main():
    # interactive
    readline.parse_and_bind("tab: complete")
    print("Insert a file name of operator: ")
    file_op = input()
    print("Insert operator rank J: ")
    op_rankJ = int(input())
    print("Insert operator rank P (1 or -1): ")
    op_rankP = int(input())
    if(op_rankP != -1 and op_rankP != 1):
        print("Parity of operator has to be 1 or -1")
        exit()
    print("Insert operator rank Z: ")
    op_rankZ = int(input())

    print("Insert a file name of transition densities: ")
    file_td = input()
    print("Insert bra state J: ")
    Jbra = int(input())
    print("Insert ket state J: ")
    Jket = int(input())
    print("Insert label of wave function for bra state: ")
    wf_label_bra = int(input())
    print("Insert label of wave function for ket state: ")
    wf_label_ket = int(input())

    print('')
    print(' main (calculation for the observable using transition density)')
    print('')

    Op = Operator.Operator(file_op, op_rankJ, op_rankP, op_rankZ)
    TD = TransitionDensity.TransitionDensity(file_td, Jbra, Jket, wf_label_bra, wf_label_ket)

    Op.read_operator_file()
    TD.read_td_file()
    TD.set_orbits(Op.orbs)
    orbs = Op.orbs

    zero = Op.zero
    one = 0.0
    for a in range(1,orbs.norbs+1):
        oa = orbs.get_orbit(a)
        for b in range(1,orbs.norbs+1):
            ob = orbs.get_orbit(b)
            if(op_rankJ == 0 and op_rankZ ==0):
                one += Op.get_obme(a,b) * TD.get_obtd(a,b,Op.rankJ,Op.rankZ) * \
                        np.sqrt(oa.j+1) / np.sqrt(2*Jbra+1)
            else:
                one += Op.get_obme(a,b) * TD.get_obtd(a,b,Op.rankJ,Op.rankZ)

    two = 0.0
    for a in range(1,orbs.norbs+1):
        for b in range(a,orbs.norbs+1):

            for c in range(1,orbs.norbs+1):
                for d in range(c,orbs.norbs+1):
                    oa = orbs.get_orbit(a)
                    ob = orbs.get_orbit(b)
                    oc = orbs.get_orbit(c)
                    od = orbs.get_orbit(d)

                    for Jab in range( int(abs(oa.j-ob.j)/2), int((oa.j+ob.j)/2)+1):
                        if(a == b and Jab%2 == 1): continue
                        for Jcd in range( int(abs(oc.j-od.j)/2), int((oc.j+od.j)/2+1)):
                            if(c == d and Jcd%2 == 1): continue
                            if(not abs(Jab-Jcd) <= Op.rankJ <= (Jab+Jcd)): continue
                            if(op_rankJ == 0 and op_rankZ ==0):
                                two += Op.get_tbme(a,b,c,d,Jab,Jcd) * TD.get_tbtd(a,b,c,d,Jab,Jcd,Op.rankJ,Op.rankZ) * \
                                        np.sqrt(2*Jab+1)/np.sqrt(2*Jbra+1)

                            else:
                                two += Op.get_tbme(a,b,c,d,Jab,Jcd) * TD.get_tbtd(a,b,c,d,Jab,Jcd,Op.rankJ,Op.rankZ)
    print('')
    print('Calculation using: ')
    print(file_op)
    print(file_td)
    print('')
    print('zero-body contribution = {:.4f}'.format(zero))
    print(' one-body contribution = {:.4f}'.format(one))
    print(' two-body contribution = {:.4f}'.format(two))
    print('                 Total = {:.4f}'.format(zero+one+two))

if(__name__ == "__main__"):
    main()
