package org.sciencework;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

public interface Permittivity {
	public double getReal( double x );
	public double getImg( double x );
	
	public DerivativeStructure getRealDerivative( DerivativeStructure x );
	public DerivativeStructure getImgDerivative( DerivativeStructure x );
}
