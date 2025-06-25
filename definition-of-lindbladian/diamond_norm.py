from lindbladian_terms import lindblad_terms
import qutip as qt
import numpy as np


J=1
h=0.5
gamma=0.1

L_dict = lindblad_terms(J=J,h=h,gamma=gamma)

L_as_list = L_dict['L_as_list']
L = L_dict['L']


d_norm_list = []


# Create output file and write header
with open('definition-of-lindbladian/output.txt', 'w') as f:
    f.write("Diamond Norm Analysis for Lindbladian Terms\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Parameters: J={J}, h={h}, gamma={gamma}\n\n")

# Define term names for identification
term_names = ['L1 (Interaction)', 'L2 (Magnetic field on qubit 1)', 
              'L3 (Magnetic field on qubit 2)', 'L4 (Decay on qubit 1)', 
              'L5 (Decay on qubit 2)']

# Calculate and write results for each term
for i, (L_i, term_name) in enumerate(zip(L_as_list, term_names)):
    # Calculate diamond norm
    d_norm = qt.dnorm(L_i,solver='SCS',verbose=True)
    
    # Append to output file
    with open('definition-of-lindbladian/output.txt', 'a') as f:
        f.write(f"Term {i+1}: {term_name}\n")
        f.write(f"Diamond Norm: {d_norm:.6f}\n")
        f.write("Matrix representation:\n")
        f.write(str(L_i.full()) + "\n")
        f.write("-" * 30 + "\n\n")
    
    # Also append to the list for later use
    d_norm_list.append(d_norm)

# Write summary at the end
with open('definition-of-lindbladian/output.txt', 'a') as f:
    f.write("Summary:\n")
    f.write("=" * 20 + "\n")
    for i, (term_name, d_norm) in enumerate(zip(term_names, d_norm_list)):
        f.write(f"Term {i+1}: {d_norm:.6f}\n")
    f.write(f"\nTotal Lindbladian Diamond Norm: {qt.dnorm(L,solver='SCS',verbose=True):.6f}\n")




    




