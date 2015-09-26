package org.sciencework;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

/**
 * Permittivity of ax^2+b+iv kind
 * @author akutlin
 *
 */
public class PermittivityImpl implements Permittivity {
	
	private final double a, b, v, bound;
	
	public PermittivityImpl( double a, double b, double v) {
		if ( a == 0 ) throw new UnsupportedOperationException("'a' should not be a zero");
		
		this.a = a;
		this.b = b;
		this.v = v;
		this.bound = (1-b)/a>0 ? Math.sqrt((1-b)/a) : 0;
	}

	public double getReal(double x) {
		if ( x > -bound && x < bound ) {
			return a*x*x + b;
		} else return 1;
	}

	public double getImg(double x) {
		if ( x > -bound && x < bound ) {
			return v;
		} else return 0;
	}

	public DerivativeStructure getRealDerivative(DerivativeStructure x) {
		double t = x.getPartialDerivative(0);
		if ( t > -bound && t < bound ) {
			return x.multiply(x).multiply(a).add(b);
		} else return x.createConstant(1);
	}

	public DerivativeStructure getImgDerivative(DerivativeStructure x) {
		return x.createConstant( getImg(x.getPartialDerivative(0)) );
	}
}
