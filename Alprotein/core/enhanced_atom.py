"""Enhanced atom class with calculation parameters."""

class EnhancedAtom:
    """Extended atom class that wraps Biopython Atom with calculation parameters.
    
    This class maintains full compatibility with Biopython Atom while adding
    calculation-specific attributes for charges and metadata.
    """
    
    def __init__(self, biopython_atom):
        """Initialize from a Biopython Atom object."""
        # Copy all attributes from biopython atom
        self.__dict__.update(biopython_atom.__dict__)
        
        # Add calculation-specific attributes
        self.q00 = None              # Ground state charge
        self.q11 = None              # First excited state charge
        self.q01 = None              # Transition charge (for couplings)
        self.is_pigment_atom = False # Flag for pigment atoms
        self.pigment_type = None     # Type of pigment this atom belongs to
        self.calculation_ready = False # Flag indicating parameters are loaded
        
        # Store reference to original atom for compatibility
        self._original_atom = biopython_atom
    
    def __getattr__(self, name):
        """Fallback to original atom for any missing attributes."""
        return getattr(self._original_atom, name)
    
    def has_calculation_parameters(self):
        """Check if atom has all required calculation parameters."""
        return (self.q00 is not None and 
                self.q11 is not None and 
                self.calculation_ready)
