##########################################################
#
# Expand_to_surface (e2s): Expands alpha-carbon-only 
# b-factor coloring to an entire surface of the protein
# of your choice.
#
# AUTHOR: Jason Vertrees -- Python code; Warren DeLano,
#         the original code.
#
# COPYRIGHT: BSDL.  Feel free to use it.
#
# DATE: 2007-12-03
#
##########################################################
from pymol import cmd
def e2s(sel):
        """
        e2s: Expand to surface

        e2s will color a surface based upon the b-factors
        from the alpha carbons in your selection.

        usage: e2s protName
        """

        cmd.create("ca_obj", sel + " and n. CA" )
        cmd.ramp_new("ramp_obj", "ca_obj", [0, 10], [-1, -1, 0] );
        cmd.set("surface_color", "ramp_obj", sel )

cmd.extend("e2s", e2s);