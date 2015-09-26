package org.sciencework;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;
import org.sciencework.core.ComplexDerivativeStructure;

public class PermittivityTest {
	
	private double getBound( double a, double b) {
		if ( a == 0 ) throw new UnsupportedOperationException("'a' should not be a zero");
		return (1-b)/a>0 ? Math.sqrt((1-b)/a) : 0;
	}
	
	@Test
	public void value() {
		double a = 2;
		double b = -5;
		double v = -0.03;
		double bound = getBound(a, b);
		
		Permittivity e = new Permittivity( a, b, v );
		for ( double i = -100; i < 100; i += 0.732 ) {
			double real = 1;
			double img = 0;
			if ( i >= -bound && i <= bound ) {
				real = a*i*i + b;
				img = v;
			}
			assertEquals( new Complex(real,img), e.value(i) );
		}
	}
	
	@Test
	public void derivative() {
		double a = 0.125;
		double b = -7;
		double v = -0.03;
		double bound = getBound(a, b);

		Permittivity e = new Permittivity( a, b, v );
		for ( double i = -100; i < 100; i += 0.232 ) {
			DerivativeStructure s = new DerivativeStructure(1,3,0,i);
			ComplexDerivativeStructure cs = new ComplexDerivativeStructure(s, s);
			
			Complex first = new Complex( e.value(cs).getReal().getPartialDerivative(1), e.value(cs).getImg().getPartialDerivative(1));
			Complex second = new Complex( e.value(cs).getReal().getPartialDerivative(2), e.value(cs).getImg().getPartialDerivative(2));
			Complex third = new Complex( e.value(cs).getReal().getPartialDerivative(3), e.value(cs).getImg().getPartialDerivative(3));
			
			Complex theFirst = new Complex(0,0);
			Complex theSecond = new Complex(0,0);
			Complex theThird = new Complex(0,0);

			
			if ( i > -bound && i < bound ) {
				theFirst = new Complex(2*a*i,0);
				theSecond = new Complex(2*a,0);
				theThird = new Complex(0,0);
			} 
//			May be sometimes...
//			else if ( i == bound || i == -bound ) {
//				theFirst = Complex.NaN;
//				theSecond = Complex.NaN;
//				theThird = Complex.NaN;
//			}
			
			assertEquals( theFirst, first);
			assertEquals( theSecond, second);
			assertEquals( theThird, third);
		}
	}

}
