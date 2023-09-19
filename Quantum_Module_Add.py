from math import pi
from projectq import MainEngine
from projectq.ops import *
from projectq.backends import CommandPrinter
from projectq.meta import Loop, Compute, Uncompute, Control

def quantumAdd(eng, a, b, n):
    QFT | b
    for i in range(0, n):
        for j in range(0, n - i):
            with Control(eng, a[n - 1 - i - j]):
                R(pi/(2**j)) | b[n - 1 - i]
    get_inverse(QFT) | b
    
def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def print_hex(eng, qubits, n):
    for i in reversed(range(n)):
        temp = 0
        temp = temp+int(qubits[4*i+3])*8
        temp = temp+int(qubits[4*i+2])*4
        temp = temp+int(qubits[4*i+1])*2
        temp = temp+int(qubits[4*i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

def print_state(eng, b, n):
    All(Measure) | b
    print('0x', end='')
    print_hex(eng, b, n)
    print('\n')

# eng = MainEngine(backend=CommandPrinter(accept_input=False))
eng = MainEngine()

size = 2
n = 4*size
a = eng.allocate_qureg(n)
b = eng.allocate_qureg(n)
c = eng.allocate_qureg(n)

Round_constant_XOR(eng, a, 0xff, n)
print_state(eng, a, size)
Round_constant_XOR(eng, b, 0xfe, n)
print_state(eng, b, size)

quantumAdd(eng, a, b, n)

print_state(eng, b, size)
