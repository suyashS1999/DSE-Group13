import numpy as np

#--------------------------------------------------------------------------------

def pe_squared(rho_inf,c,P,D,F,r,M,theta):
    pe = np.zeros((len(F),len(theta)))
    for i in range(len(F)):
        pe[i,:] = (rho_inf*c*P*D*F[i])/(4*np.pi*(r**2)*(1-M*np.cos(np.radians(theta)))**4)

    return pe

#--------------------------------------------------------------------------------

def powers(K,M,a,G,rho_inf,c,bw):

    p = K*((M)**a)*G*(rho_inf*c**3 *bw**2)

    return p

#--------------------------------------------------------------------------------

def Strouhal_number(f,L,M,theta,c):
    S = np.zeros((len(f),len(theta)))
    for i in range(len(f)):
        S[i,:] = (f[i]*L*(1-M*np.cos(np.radians(theta))))/(M*c)

    return S

#--------------------------------------------------------------------------------

def geometry_func(A,b,rho_inf,M,c,dynamic_visc,df,number_of_wheels,d,noise_source):

    if noise_source == "Clean wing":

        G = 0.37*(A/(b**2))*((rho_inf*M*c*A)/(dynamic_visc*b))**(-0.2)

    elif noise_source == "Slats":

        G = 0.37*(A/(b**2))*((rho_inf*M*c*A)/(dynamic_visc*b))**(-0.2)

    elif noise_source == "Flaps":

        G = (A/(b**2))*np.sin(np.radians(df))**2

    elif noise_source == "Main landing gear":

        G = number_of_wheels*(d/b)**2

    else:

        print("Incorrect noise source (L)")

        G = 0

    return G

#--------------------------------------------------------------------------------

def length_scale(G,b,A,d,noise_source):

    if noise_source == "Clean wing":

        L = G*b

    elif noise_source == "Slats":

        L = G*b

    elif noise_source == "Flaps":

        L = A/b

    elif noise_source == "Main landing gear":

        L = d

    else:

        print("Incorrect noise source (L)")

        L = 0

    return L

#--------------------------------------------------------------------------------

def constants_Ka(noise_source,number_of_wheels):

    if noise_source == "Clean wing":

        K = 4.464E-5

        a = 5

    elif noise_source == "Slats":

        K = 4.464E-5

        a = 5

    elif noise_source == "Flaps":

        K = 2.787E-4

        a = 6

    elif noise_source == "Main landing gear":

        if number_of_wheels <= 2:

            K = 4.349E-4

        elif number_of_wheels == 4:

            K = 3.414E-4

        else:

            print("Method does not account for this amount of wheels")

            K = 3.414E-4

        a = 6

    else:

        print("Incorrect noise source (K,a)")

        K = 0

        a = 0
    
    return K,a

#--------------------------------------------------------------------------------

def directivity(theta,phi,df,noise_source):

    if noise_source == "Clean wing":

        D = 4*(np.cos(np.radians(phi))**2)*(np.cos(np.radians(theta/2))**2)

    elif noise_source == "Slats":

        D = 4*(np.cos(np.radians(phi))**2)*(np.cos(np.radians(theta/2))**2)

    elif noise_source == "Flaps":

        D = 3*(np.sin(np.radians(df))*np.cos(np.radians(theta)) + np.cos(np.radians(df))*np.sin(np.radians(theta))*np.cos(np.radians(phi)))**2

    elif noise_source == "Main landing gear":

        D = (3/2)*np.sin(np.radians(theta))**2

    else:

        print("Incorrect noise source (D)")

        D = 0
    
    return D

#--------------------------------------------------------------------------------

def dimensionless_empirical_spectral_func(S, number_of_wheels, noise_source):

    if noise_source == "Clean wing":

        F = (0.613*(10*S)**4)*(((10*S)**1.5 + 0.5)**(-4))

    elif noise_source == "Slats":

        F = (0.613*(10*S)**4)*(((10*S)**1.5 + 0.5)**(-4)) + (0.613*(2.19*S)**4)*(((2.19*S)**1.5 + 0.5)**(-4))

    elif noise_source == "Flaps":
        F = np.zeros(len(S))
        for i in range(len(S)):
            
            if S[i] < 2:

                F[i] = 0.048*S[i]

            elif 2 <= S[i] <= 20:
                    
                F[i] = 0.1406*((S[i])**(-0.55))
                

            elif S[i] > 20:

                F = 216.49*((S[i])**(-3))

            else:

                print("Something is wrong in F(S)")

                F[i,j] = 0

    elif noise_source == "Main landing gear":

        if number_of_wheels <= 2:

            F = 13.59*(S**2)*(S**2 + 12.5)**(-2.25)

        elif number_of_wheels == 4:

            F = 0.0577*(S**2)*(0.25*(S**2) + 1)**(-1.5)

        else:

            print("Method does not account for this amount of wheels")

            F = 0

    else:

        print("Incorrect noise source F(S)")

        F = 0

    return F

#--------------------------------------------------------------------------------

def SPL(pe_sq,pe0):

    SPL = 10*np.log10((pe_sq)/(pe0**2))

    return SPL

 #--------------------------------------------------------------------------------

def dL_A(f):

    dL_A = -145.528 + 98.262*np.log10(f) - 19.509*(np.log10(f)**2) + 0.975*(np.log10(f)**3)

    return dL_A

 #--------------------------------------------------------------------------------

def L_A(SPL,dL_A):
    L_A = np.zeros((len(dL_A),len(SPL)))

    for i in range(len(SPL)):
    
        L_A[:,i] = 10*np.log10(np.sum(10**((SPL[:,i]+dL_A)/10)))

    return L_A

 #--------------------------------------------------------------------------------