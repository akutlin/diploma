package org.sciencework;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

public class Calculator {

	public static void main(String[] args) {
		
		double a = 1;
		double b = -1;
		double v = -0.1;
		
		PermittivityImpl e = new PermittivityImpl( a, b, v );
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1,0);
		double teta = 0;
		
		TE wave = new TE(teta, e, 1E11);
		
		AdamsMoultonIntegrator integrator = new  AdamsMoultonIntegrator(3, 0.0000001, 0.001, 0.001, 0.001);
		
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
		
		double dissipation = ( S1N - S0N )/S1N; //S1 is the begining!!!
		
		System.out.println( (double) Math.round( dissipation * 1000 ) / 1000 );
		
	}

}
