package org.sciencework;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.nonstiff.AdamsBashforthIntegrator;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
//import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;

public class Calculator {

	public static void main(String[] args) {
		
		double a = 1;
		double b = -1;
		double v = -0.1;
		
		Permittivity e = new Permittivity( a, b, v );
		
		Complex Ey0 = new Complex(1,0);
		Complex Bz0 = new Complex(1,0);
		double teta = 0;
		
		TE wave = new TE();
		wave.setTeta(teta);
		wave.setPermittivity(e);
		wave.setFrequency(100000000000L);
		
		AdamsMoultonIntegrator integrator = new  AdamsMoultonIntegrator(3, 0.0000001, 0.001, 0.001, 0.001);
		
		double y[] = new double[] { Ey0.getReal(), Ey0.getImaginary(), Bz0.getReal(), Bz0.getImaginary() };
		
		double t0 = 2;
		double t1 = -2;
		integrator.integrate(wave, t0, y, t1, y);
		
		Complex Ey1 = new Complex( y[0], y[1] );
		Complex Bz1 = new Complex( y[2], y[3] );
		
		double Sx0 = Ey0.multiply(Bz0.conjugate()).getReal();
		double Sz0 = Ey0.multiply( Bz0.multiply( Math.tan(teta)).conjugate() ).getReal();
		double S0 = Math.sqrt( Math.pow(Sx0, 2) + Math.pow(Sz0, 2) );
		
		double Sx1 = Ey1.multiply(Bz1.conjugate()).getReal();
		double Sz1 = Ey1.multiply( Bz1.multiply( Math.tan(teta)).conjugate() ).getReal();
		double S1 = Math.sqrt( Math.pow(Sx1, 2) + Math.pow(Sz1, 2) );
		
		double dissipation = ( S1 - S0 )/S1; //S1 is the begining!!!
		
		System.out.println( (double) Math.round( dissipation * 1000 ) / 1000 );
		
	}

}
