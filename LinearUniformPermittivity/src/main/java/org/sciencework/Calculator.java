package org.sciencework;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

public class Calculator {

	public static void main(String[] args) {
		
		double a = 1;
		double b = -1;
		double v = -0.1;
		double teta = 0;
		double k = 100;
		
		PermittivityImpl e = new PermittivityImpl( a, b, v );
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1 * Math.cos(teta),0);
		
		TE wave = new TE(teta, e, k * TE.c);
		
		DormandPrince853Integrator integrator = new  DormandPrince853Integrator(1E-11, 1E-3, 1E-10, 1E-10);
		
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = 2;
		double t1 = -2;
		integrator.integrate(wave, t0, y, t1, y);
		
		Complex Ey1 = new Complex( y[0], y[1] );
		Complex Bz1 = new Complex( y[2], y[3] );
		
		RealVector S0 = wave.getPowerFlow(Ey0, Bz0);
		double S0N = S0.getNorm();
		
		RealVector S1 = wave.getPowerFlow(Ey1, Bz1);
		double S1N = S1.getNorm();
		
		double dissipation = ( S0N - S1N )/S1N; //S1 is the begining!!!
		
		System.out.println( (double) Math.round( dissipation * 100 ) / 100 );
		
	}

}
