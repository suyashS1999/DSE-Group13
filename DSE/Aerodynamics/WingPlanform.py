from numpy import*

#%% ------------------- Input data -------------------
# Global Main Wing
AR = 17;									# Wing Aspect ratio [-]
S = 142.520611;								# Wing surface area [m^2]
span = 49.22245815;							# Wing span [m]
ave_chord = S/span;							# Average chord lenght [m]
taper_ratio = 0.3;							# Taper ratio [-]
# Inboard Wing
span_inboard = 36;							# Span of inboard wing [m]
span_outboard = span - span_inboard;		# Span of outboard wing [m]

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

root_chord, tip_chord = Calc_root_tip_chordMain(span, AR, taper_ratio);
inb, outb = InboardOutboard_wing_parms(S, span_inboard, span_outboard, root_chord, tip_chord);
print(inb, outb);