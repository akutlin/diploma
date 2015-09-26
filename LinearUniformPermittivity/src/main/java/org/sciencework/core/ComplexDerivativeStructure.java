package org.sciencework.core;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

/**
 * It is a complex number adapter for {@link DerivativeStructure}
 * @author akutlin
 *
 */

public class ComplexDerivativeStructure {
	
	private final DerivativeStructure real, img;
	
	public ComplexDerivativeStructure( DerivativeStructure real, DerivativeStructure img ) {
		this.real = real;
		this.img = img;
	}
	
	public DerivativeStructure getReal() {
		return real;
	}
	
	public DerivativeStructure getImg() {
		return img;
	}
}
