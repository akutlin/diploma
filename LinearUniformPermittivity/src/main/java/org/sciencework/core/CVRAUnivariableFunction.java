package org.sciencework.core;

import org.apache.commons.math3.complex.Complex;

/**
 * CVRA means Complex Value and Real Argument
 * It is a complex number adapter for {@link UnivariableFunction}
 * @author akutlin
 * 
 */
public interface CVRAUnivariableFunction {
	public Complex value( double arg );
}
