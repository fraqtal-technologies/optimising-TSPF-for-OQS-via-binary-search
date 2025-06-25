import qutip as qt
import numpy as np

def lindblad_terms(J,h,gamma):
    
    """
    Compute Lindbladian terms for a two-qubit system with interaction, magnetic field, and decay.
    
    Parameters:
    -----------
    J : float
        Strength of interaction/coupling between qubits (real number)
    h : float  
        Magnetic field parameter (real number)
    gamma : float
        Decay rate (real number)
        
    Returns:
    --------
    dict
        Dictionary containing the total Lindbladian operator 'L' and 
        list of individual terms 'L_as_list'
    """
    
    
    L1 = qt.liouvillian(H=1j*J*qt.tensor(qt.sigmaz(),qt.sigmaz()),
                        c_ops=[])
    
    L2 = qt.liouvillian(H=1j*h*qt.tensor(qt.sigmax(),qt.qeye(2)),
                        c_ops=[])
    
    L3 = qt.liouvillian(H=1j*h*qt.tensor(qt.qeye(2),qt.sigmax()),
                        c_ops=[])
    
    L4 = gamma*qt.liouvillian(H=0*qt.tensor(qt.qeye(2),qt.qeye(2)),
                        c_ops=[qt.tensor(qt.sigmaz(),qt.qeye(2))])
    
    L5 = gamma*qt.liouvillian(H=0*qt.tensor(qt.qeye(2),qt.qeye(2)),
                        c_ops=[qt.tensor(qt.qeye(2),qt.sigmaz())])
    
    return {"L":L1+L2+L3+L4+L5,
            "L_as_list":[L1,L2,L3,L4,L5]
            }
    




