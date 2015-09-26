package org.sciencework;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import javax.annotation.Nonnull;

public class TE implements FirstOrderDifferentialEquations {
	
	private double sinOfTetaSquare;
	private double k;
	private Permittivity e;
	
	private final int DIM = 4;
	public final long c = 10000000000L;
	
	public void setTeta( double teta ) {
		sinOfTetaSquare = Math.pow( Math.sin( teta ), 2);
	}
	
	public void setPermittivity( Permittivity e ) {
		this.e = e;
	}
	
	public void setFrequency( double w ) {
		k = w / c;
	}

	public int getDimension() {
		return DIM;
	}

	/**
	 * y = { re(Ey), im(Ey), re(Bz), im(Bz) }
	 */
	public void computeDerivatives(double t, @Nonnull double[] y, @Nonnull double[] yDot)
			throws MaxCountExceededException, DimensionMismatchException {
		if ( y.length != DIM ) throw new DimensionMismatchException( y.length, DIM);
		if ( yDot.length != DIM ) throw new DimensionMismatchException( yDot.length, DIM);
		
		Complex E = e.value(t);
		double Er = E.getReal();
		double Ei = E.getImaginary();
		
		Complex Ey1 = new Complex( y[0], y[1] );
		Complex Bz1 = new Complex( y[2], y[3] );
		
		double Sx0 = Ey1.multiply(Bz1.conjugate()).getReal();
		double Sz0 = Ey1.multiply( Bz1.multiply( Math.tan( Math.asin( Math.pow(sinOfTetaSquare, 2)))).conjugate() ).getReal();
		double S0 = Math.sqrt( Math.pow(Sx0, 2) + Math.pow(Sz0, 2) );
//		System.out.println( "(" + y[0] + "," + y[1] + "), (" + y[2] + "," + y[3] + ")" );
		System.out.println( E + "-> ( " + Sx0 + ", " + Sz0 + "): " + S0);
//		System.out.println(S0);
				
		yDot[0] = -k * y[3];
		yDot[1] = k * y[2];
		yDot[2] = -k * ( Ei * y[0] + ( Er - sinOfTetaSquare ) * y[1] );
		yDot[3] = k * ( y[0] * ( Er - sinOfTetaSquare ) - Ei * y[1] );
	}

}
