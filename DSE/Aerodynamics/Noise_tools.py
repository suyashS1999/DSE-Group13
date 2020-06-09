import numpy as np

#--------------------------------------------------------------------------------

def pe_squared(rho_inf,c,P,D,F,r,M,theta):

    pe = (rho_inf*c*p*D*F)/(4*np.pi()*(r**2)*(1-M*np.cos(theta))**4)

    return pe

#--------------------------------------------------------------------------------

def power(K,M,a,G):

    p = K*(M**a)*G

    return pgit

#--------------------------------------------------------------------------------

def Strouhal_number(f,L,M,theta,c):

    S = (f*L*(1-M*np.cos(theta)))/(M*c)

    return S

#--------------------------------------------------------------------------------

def geometry_func(A,b,rho_inf,M,c,dynamic_visc,df,number_of_wheels,d,noise_source):

    if noise_source == "Clean wing":

        G = 0.37*(A/(b**2))*((rho_inf*M*c*A)/(dynamic_visc*b))**(-0.2)

    elif noise_source == "Slats":

        G = 0.37*(A/(b**2))*((rho_inf*M*c*A)/(dynamic_visc*b))**(-0.2)

    elif noise_source == "Flaps":

        G = (A/(b**2))*np.sin(df)**2

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

        D = 4*(np.cos(phi)**2)*(np.cos(theta/2)**2)

    elif noise_source == "Slats":

        D = 4*(np.cos(phi)**2)*(np.cos(theta/2)**2)

    elif noise_source == "Flaps":

        D = 3*(np.sin(df)*np.cos(theta) + np.cos(df)*np.sin(theta)*np.cos(phi))**2

    elif noise_source == "Main landing gear":

        D = (3/2)*np.sin(theta)**2

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

        if S < 2:

            F = 0.048*S

        elif 2 <= S <= 20:

            F = 0.1406*(S**(-0.55))

        elif S > 20:

            F = 216.49*(S**(-3))

        else:

            print("Something is wrong in F(S)")

            F = 0

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