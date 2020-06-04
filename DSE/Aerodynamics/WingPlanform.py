from numpy import*
#%% ------------------- Functions -------------------
def Calc_root_tip_chordMain(span, AR, taper_ratio):
	""" Function to compute the root and tip chord for the main wing
	Input:
		span = Span of the main wing [m]
		AR = The value for the AR of the main wing [m]
		taper_ratio = Taper ratio of the main wing [-]
	Output:
		root_chord = Root Chord [m]
		tip_chord = Tip Chord [m]
	"""
	half_span = span/2;		half_AR = AR/2;
	root_chord = 2*half_span/(half_AR*(1 + taper_ratio));
	tip_chord = taper_ratio*root_chord;
	return root_chord, tip_chord;

def InboardOutboard_wing_parms(S, span_inboard, span_outboard, root_chord, tip_chord):
	""" Function to compute the wing parameters for the inboard and outboard wings
	Input:
		S = Main Wing area [m^2]
		span_inboard = Span of inboard wing [m]
		span_outboard = Span of outboard wing [m]
		root_chord = Root chord of main wing [m]
		tip_chord = tip chord of main wing [m]
	Output:
		(root_chord_inb, tip_chord_inb, S_inb, AR_inb, taper_ratio_inb) = Inboard wing parameters [SI units] (tuple)
		(root_chord_outb, tip_chord_outb, S_outb, AR_outb, taper_ratio_outb) = Outboard wing parameters [Si units] (tuple)
	"""
	root_chord_inb = root_chord;
	tip_chord_inb = (S/2 - span_inboard*root_chord/4 - span_outboard*tip_chord/4)/((span_inboard + span_outboard)/4);
	S_inb = span_inboard/2*(root_chord_inb + tip_chord_inb);
	AR_inb = span_inboard**2/S_inb;
	taper_ratio_inb = tip_chord_inb/root_chord_inb;

	root_chord_outb = tip_chord_inb;
	tip_chord_outb = tip_chord;
	S_outb = span_outboard/2*(root_chord_outb + tip_chord_outb);
	AR_outb = span_outboard**2/S_outb;
	taper_ratio_outb = tip_chord_outb/root_chord_outb;

	return (root_chord_inb, tip_chord_inb, S_inb, AR_inb, taper_ratio_inb), (root_chord_outb, tip_chord_outb, S_outb, AR_outb, taper_ratio_outb);

def Compute_MAC(root_chord, tip_chord, sweep_LE, span):
	""" Function to compute the Mean Aerodynamic chord
	Input:
		root_chord = Root chord of the main wing [m]
		tip_chord = Tip chord of main wing [m]
		sweep_LE = leading edge sweep [rad]
		span = wing span [m]
	Output:
		MAC = Mean Aerodynamic Chord [m]
		MAC_y = span location of MAC [m]
	"""
	a = lambda x: x*tan(sweep_LE);
	f = lambda x: tip_chord + root_chord/2 + ((a(span/2) - tip_chord/2 - root_chord/2)/(span/2))*x;
	x = (tip_chord + root_chord/2)/((1.5*tip_chord + 1.5*root_chord)/(span/2));
	a_prime = a(x);
	c = f(x);
	MAC = (c - a_prime - tip_chord)*2;
	MAC_y = x;
	return MAC, MAC_y;


def Optimal_taper(sweep_quart):
	""" Function to compute the optimum taper ratio
	Input:
		sweep_quart = Quater chord sweep angle [rad]
	Output:
		taper = optimum taper ratio
	"""
	taper = 0.45*exp(-0.036*sweep_quart);
	return taper