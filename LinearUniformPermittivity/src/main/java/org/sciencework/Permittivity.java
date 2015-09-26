package org.sciencework;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.sciencework.core.CVRAUnivariateDifferentiableFunction;
import org.sciencework.core.ComplexDerivativeStructure;

/**
 * Permittivity of ax^2+b+iv kind
 * @author akutlin
 *
 */
public class Permittivity implements CVRAUnivariateDifferentiableFunction {
	
	private final double a, b, v, bound;
	private final Real realPart = new Real();
	private final Img imgPart = new Img();
	
	private class Real implements UnivariateDifferentiableFunction {
		
		public double value(double x) {
			if ( x > -bound && x < bound ) {
				return a*x*x + b;
			} else return 1;
		}

		public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {
			double x = t.getPartialDerivative(0);
			if ( x > -bound && x < bound ) {
				return t.multiply(t).multiply(a).add(b);
			} else return t.createConstant(1);
		}
		
	}
	
	private class Img implements UnivariateDifferentiableFunction {
		
		public double value(double x) {
			if ( x > -bound && x < bound ) {
				return v;
			} else return 0;
		}

		public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {
			return t.createConstant( value(t.getPartialDerivative(0)) );
		}
		
	} 
	
	public Permittivity( double a, double b, double v) {
		if ( a == 0 ) throw new UnsupportedOperationException("'a' should not be a zero");
		
		this.a = a;
		this.b = b;
		this.v = v;
		this.bound = (1-b)/a>0 ? Math.sqrt((1-b)/a) : 0;
	}

	public Complex value(double arg) {
		double real = realPart.value(arg);
		double img = imgPart.value(arg);
		return new Complex(real,img);
	}

	public ComplexDerivativeStructure value(ComplexDerivativeStructure t) throws DimensionMismatchException {
		DerivativeStructure real = realPart.value( t.getReal() );
		DerivativeStructure img = imgPart.value( t.getImg() );
		return new ComplexDerivativeStructure(real, img);
	}
}
