import numpy as np
import matplotlib.pyplot as plt

pi = np.pi


# Lift coeficcient
def LiftValue(alpha):
	if alpha<=0.25:
		CL = 2*np.pi*alpha
	else:
		CL = np.pi/2*np.cos(alpha)/np.cos(0.25)
	return CL

# Drag coeficcient
def DragValue(alpha):
	if alpha<=0.25:
		CD = 0.224*alpha**2 + 0.006
	elif alpha>=0.25 and alpha<= 0.30  :
		CD = 16.6944*alpha**2 - 1.0234
	else:
		CD = np.pi/2*np.sin(alpha)/np.sin(0.25)
	return CD

# Array of force coefficients for an array of angles of attack
def ArrayForceCoefficients(alpha):

	CL = []; CD = []
	for i in alpha:
		CL.append(LiftValue(i))
		CD.append(DragValue(i))

	CL = np.asarray(CL)
	CD = np.asarray(CD)
	
	return CL, CD



# Derivation for calculating the induced angle
def depsiloni(Cbb,rR,CL,k,betaTip,epsiloni,epsiloninf):

	Term1 = Cbb/(8*rR)*CL

	Term2 = -np.arccos(np.exp(-k*(1-rR)/(2*np.sin(betaTip))))

	Term3 = np.tan(epsiloni)*np.sin(epsiloni+epsiloninf)

	Term4 = np.tan(epsiloni)*np.cos(epsiloni+epsiloninf) \
			+ 1/np.cos(epsiloni)**2 * np.sin(epsiloni+epsiloninf)

	fx = Term1 + Term2*Term3
	dfx = Term2*Term4

	return fx, dfx


zeroliftangle=-2.1
zeroliftangle=zeroliftangle*np.pi/180

# Radius
Nn = 30
rR = np.linspace(0.1,0.99,Nn)
print(rR)

cbdp = 0.075*np.sqrt(1-rR**2)


# Aerodynamic pitch-to-diameter ratio
Kc = 0.5
# Local pitch-to-diameter ratio
K = pi*rR* ( Kc - pi*rR*np.tan(zeroliftangle) ) / ( pi*rR + Kc*np.tan(zeroliftangle) )

# Aerodynamic pitch angle 
beta = np.arctan(K/(pi*rR))
betaTip = beta[-1] #Last element of the arrow

# print(beta)
# print(beta*180/pi)

# Advance angle 
J = 0.25
epsiloninf = np.arctan(J/(pi*rR) )

epsiloni = np.linspace(1,1,30)
# CLL = np.linspace(0,0,30)

# Number of blades
k = 2
# Chord length ratio
Cbb = k * 0.075*np.sqrt(1-rR**2)

# lamb = 2*np.pi*r * (lambc - 2*np.pi*r*tan(zeroliftangle)) / (2*np.pi*r + lambc*tan(zeroliftangle))


alphab = beta - epsiloninf - epsiloni

CL, CD = ArrayForceCoefficients(alphab)


print('++++++++++++++++++')

epsiloni0 = CL*0.01

for i in  range(50):

	fx, dfx = depsiloni(Cbb,rR,CL,k,betaTip,epsiloni0,epsiloninf)
	epsiloni1 = epsiloni0 - fx/dfx
	print(i)

	alphab = beta - epsiloninf - epsiloni1

	# Calculate new force coefficients
	CL, CD = ArrayForceCoefficients(alphab)

	Error = abs(sum(epsiloni1-epsiloni0))
	print(Error)
	epsiloni0 = epsiloni1

	if Error < 10e-6: break


print(epsiloni1*180/pi)


# print(1/np.cos(epsiloni))

# sec^2(x) sin(x + y) + cos(x + y) tan(x)

# plt.plot(rR, cbdp, 'bo')
plt.plot(rR, epsiloni1*180/pi, 'ro')
plt.show()