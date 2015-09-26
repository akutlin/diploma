package org.sciencework.core;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

public class CVRAFirstOrderIntegratorImpl implements CVRAFirstOrderIntegrator {

	private FirstOrderIntegrator integrator;
	
	public double integrate(CVRAFirstOrderDifferentialEquations equations, double t0, Complex[] y0, double t, Complex[] y) {
		
		FirstOrderDifferentialEquations real = equations.getRealEquations();
		FirstOrderDifferentialEquations img = equations.getImaginaryEquations();
		
		int dim = equations.getDimension();
		double[] realY0 = new double[ dim ];
		double[] imgY0 = new double[ dim ];
		double[] realY = new double[ dim ];
		double[] imgY = new double[ dim ];
		
		for ( int i = 0; i < dim; i++ ) {
			realY0[i] = y0[i].getReal();
			imgY0[i] = y0[i].getImaginary();
		}
		
		double timeR = integrator.integrate(real, t0, realY0, t, realY);
		double timeI = integrator.integrate(img, t0, imgY0, t, imgY);
		
		for( int i = 0; i < dim; i++ ) {
			y[i] = new Complex( realY[i], imgY[i] );
		}

		return Math.min( timeR, timeI );
	}

	public void setFirstOrderIntegrator(FirstOrderIntegrator integrator) {
		this.integrator = integrator;
	}

}
