package org.sciencework;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import javax.annotation.Nonnull;

public class TE implements FirstOrderDifferentialEquations {
	
	private final double sinOfTetaSquare;
	private final double k;
	private final Permittivity e;
	
	private final int DIM = 4;
	public static final long c = (long) 1E10;
	
	public TE( double teta, Permittivity e, double w ) {
		sinOfTetaSquare = Math.pow( Math.sin( teta ), 2);
		this.e = e;
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
		
		double Er = e.getReal(t);
		double Ei = e.getImg(t);
		
//		RealVector S = getPowerFlow( new Complex( y[0], y[1] ), new Complex( y[2], y[3] ));
//		double S0 = S.getNorm();
//		System.out.println( new Complex(Er,Ei) + " -> " + S + " : " + S0);
//		System.out.println(S0);
//		System.out.println( t + " -> (" + y[0] + "," + y[1] + "), (" + y[2] + "," + y[3] + ") : " + S0 );
				
		yDot[0] = -k * y[3];
		yDot[1] = k * y[2];
		yDot[2] = -k * ( Ei * y[0] + ( Er - sinOfTetaSquare ) * y[1] );
		yDot[3] = k * ( y[0] * ( Er - sinOfTetaSquare ) - Ei * y[1] );
	}
	
	/**
	 * Pointing's vector
	 * @param Ey
	 * @param Bz
	 * @return { Sx, Sz }
	 */
	public RealVector getPowerFlow( Complex Ey, Complex Bz ) {
		double Sx = Ey.multiply( Bz.conjugate()).getReal();
		double Sz = Ey.multiply( Bz.multiply( Math.tan( Math.asin( Math.sqrt(sinOfTetaSquare)))).conjugate() ).getReal();
		return MatrixUtils.createRealVector( new double[] { Sx, Sz} );
	}

}
